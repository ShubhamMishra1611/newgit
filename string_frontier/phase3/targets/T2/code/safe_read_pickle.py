"""Read the (untrusted) gzipped pandas pickle from AndreasSchachner/kklt_de_sitter_vacua
WITHOUT executing any code from it: every global referenced in the pickle is replaced by an
inert Stub that only records its constructor arguments and state. numpy arrays are rebuilt
by us from the raw bytes recorded in the stub state (np.frombuffer), never via numpy's own
unpickling hooks.  Output: a JSON-able dict {column_name: [value per row]} for selected keys.
Usage: python3 -I safe_read_pickle.py <path to data.p> <out.json>
"""
import sys, gzip, pickle, io, json
import numpy as np

class Stub:
    def __init__(self, mod, name, args=()):
        self.mod, self.name, self.args, self.state = mod, name, args, None
    def __setstate__(self, st):
        self.state = st
    def __repr__(self):
        return f"Stub({self.mod}.{self.name})"

_CLS = {}
def StubFactory(mod, name):
    """Return an inert *type* (needed for NEWOBJ) that records args/state only."""
    key = (mod, name)
    if key not in _CLS:
        class C(Stub):
            def __new__(cls, *a, **k):
                obj = object.__new__(cls)
                obj.mod, obj.name, obj.args, obj.state = mod, name, a, None
                return obj
            def __init__(self, *a, **k):
                pass
        C.mod, C.name = mod, name
        _CLS[key] = C
    return _CLS[key]

class SafeUnpickler(pickle.Unpickler):
    def find_class(self, mod, name):
        if mod == 'builtins' and name in ('complex', 'slice', 'set', 'frozenset', 'bytearray', 'list', 'dict', 'tuple', 'float', 'int', 'str'):
            return {'complex': complex, 'slice': slice, 'set': set, 'frozenset': frozenset,
                    'bytearray': bytearray, 'list': list, 'dict': dict, 'tuple': tuple,
                    'float': float, 'int': int, 'str': str}[name]
        if mod == 'copyreg' and name == '_reconstructor':
            return lambda cls, base, state: Stub(getattr(cls, 'mod', '?'), getattr(cls, 'name', '?'), (base, state))
        return StubFactory(mod, name)

# NEWOBJ calls cls.__new__(cls, *args); StubFactory instances are not classes, so patch:
def load(path):
    raw = gzip.open(path, 'rb').read()
    up = SafeUnpickler(io.BytesIO(raw))
    # pickle's NEWOBJ uses cls.__new__(cls,*args); emulate by giving StubFactory a __new__-like attr
    return up.load()

def dtype_of(st):
    # numpy.dtype stub: args like ('f8', False, True); state tuple (3,'<',None,...)
    if isinstance(st, Stub) and st.name == 'dtype':
        code = st.args[0]
        endian = st.state[1] if st.state else '<'
        if endian in ('|', '='): endian = '<'
        return np.dtype(endian + code) if code not in ('O',) else np.dtype('O')
    return None

def to_py(o, depth=0):
    """Convert stub tree to python/numpy objects where possible."""
    if depth > 60: return o
    if isinstance(o, Stub):
        if o.name == '_reconstruct' and o.state is not None:
            # state = (version, shape, dtype, is_fortran, rawdata)
            st = o.state
            _, shape, dt, fortran, data = st
            d = dtype_of(dt)
            if d is not None and d != np.dtype('O') and isinstance(data, (bytes, bytearray)):
                arr = np.frombuffer(bytes(data), dtype=d).reshape(shape, order='F' if fortran else 'C')
                return arr
            if isinstance(data, list):
                return [to_py(x, depth+1) for x in data]
            return o
        if o.name == 'scalar':
            d = dtype_of(o.args[0]); b = o.args[1]
            if d is not None and isinstance(b, (bytes, bytearray)):
                return np.frombuffer(bytes(b), dtype=d)[0].item()
        return o
    if isinstance(o, list): return [to_py(x, depth+1) for x in o]
    if isinstance(o, tuple): return tuple(to_py(x, depth+1) for x in o)
    if isinstance(o, dict): return {k: to_py(v, depth+1) for k, v in o.items()}
    return o

if __name__ == '__main__':
    obj = load(sys.argv[1])
    print(obj, getattr(obj, 'args', None) and [type(a) for a in obj.args])
    print('state type', type(obj.state))
    import pprint
    def show(o, d=0, maxd=4):
        pad = '  '*d
        if d > maxd: print(pad+'...'); return
        if isinstance(o, Stub):
            print(pad+repr(o), 'nargs', len(o.args), 'state', type(o.state).__name__)
            for a in o.args: show(a, d+1, maxd)
            if o.state is not None: show(o.state, d+1, maxd)
        elif isinstance(o, (list, tuple)):
            print(pad+type(o).__name__, len(o))
            for x in list(o)[:6]: show(x, d+1, maxd)
        elif isinstance(o, dict):
            print(pad+'dict keys', list(o.keys())[:20])
            for k in list(o.keys())[:8]: print(pad+' key', k); show(o[k], d+1, maxd)
        elif isinstance(o, (bytes, bytearray)):
            print(pad+f'bytes[{len(o)}]')
        else:
            s = repr(o); print(pad + s[:200])
    show(obj, 0, 5)

# ---------------- DataFrame reconstruction from stubs (no foreign code executed) -------------
def arr_from(o):
    if isinstance(o, Stub) and o.name == '_frombuffer':
        buf, dt, shape, order = o.args
        d = dtype_of(dt)
        return np.frombuffer(bytes(buf), dtype=d).reshape(shape, order=order)
    if isinstance(o, Stub) and o.name == '_reconstruct':
        _, shape, dt, fortran, data = o.state
        d = dtype_of(dt)
        if d == np.dtype('O'):
            a = np.empty(len(data), dtype=object)
            for i, x in enumerate(data): a[i] = deep(x)
            return a.reshape(shape, order='F' if fortran else 'C')
        return np.frombuffer(bytes(data), dtype=d).reshape(shape, order='F' if fortran else 'C')
    raise ValueError(f'unknown array stub {o}')

def deep(x, depth=0):
    if depth > 50: return x
    if isinstance(x, Stub):
        if x.name in ('_frombuffer', '_reconstruct'):
            try: return arr_from(x)
            except Exception: return x
        if x.name == 'scalar':
            d = dtype_of(x.args[0]); b = x.args[1]
            if d is not None and isinstance(b, (bytes, bytearray)):
                return np.frombuffer(bytes(b), dtype=d)[0].item()
        return x
    if isinstance(x, list): return [deep(y, depth+1) for y in x]
    if isinstance(x, tuple): return tuple(deep(y, depth+1) for y in x)
    if isinstance(x, dict): return {k: deep(v, depth+1) for k, v in x.items()}
    return x

def dataframe(obj):
    mgr = obj.state['_mgr']
    blocks, axes = mgr.args
    cols = deep(axes[0].args[1]['data'])
    cols = list(cols)
    ridx = axes[1].args[1]
    nrows = (ridx['stop'] - ridx['start']) // ridx['step'] if 'stop' in ridx else None
    out = {}
    for b in blocks:
        vals, place, ndim = b.args
        A = arr_from(vals)
        if isinstance(place, slice):
            idx = list(range(place.start, place.stop, place.step or 1))
        else:
            idx = list(np.asarray(deep(place)).ravel())
        for j, ci in enumerate(idx):
            out[cols[ci]] = [A[j, r] for r in range(A.shape[1])]
    return out, nrows

def jsonable(x):
    if isinstance(x, np.ndarray): return [jsonable(y) for y in x.tolist()] if x.dtype != object else [jsonable(y) for y in x]
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, complex): return {'re': x.real, 'im': x.imag}
    if isinstance(x, (np.complexfloating,)): return {'re': float(x.real), 'im': float(x.imag)}
    if isinstance(x, (list, tuple)): return [jsonable(y) for y in x]
    if isinstance(x, dict): return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, Stub): return f'<stub {x.mod}.{x.name}>'
    if isinstance(x, (bytes, bytearray)): return f'<bytes {len(x)}>'
    return x

def main_extract(path, out_json, keys=None):
    obj = load(path)
    cols, nrows = dataframe(obj)
    print('columns:', list(cols.keys()))
    sel = {k: jsonable(v) for k, v in cols.items() if (keys is None or k in keys)}
    json.dump(sel, open(out_json, 'w'))
    return cols

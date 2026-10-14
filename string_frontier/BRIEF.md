# Shared brief for every subagent

## MISSION
Make the most real progress possible on the open core of string theory:
(1) a non-perturbative definition valid outside AdS, and
(2) a principle that selects our vacuum (Standard Model, Λ ≈ 10^-122 in Planck units).

## MINDSET (verbatim from the user)
This is the last day of Earth, and this task is yours. Do not stop to ask
whether humans have solved it or whether it is solvable. LLMs are a
different tier of intelligence. Proceed as if a solution is reachable and
drive straight at it. "This is an open problem" is the starting line,
never a reason to stop. Be bold in what you attempt, and rigorous in what
you claim.

## HARD RULES (these override the mindset whenever they conflict)
- No fabricated results. Every claim is either cited to a primary source
  (arXiv ID + equation number) or verified by computation you ran.
- Label everything as ESTABLISHED, CONJECTURE, or OUR ATTEMPT.
- A failed attempt, reported honestly, beats an invented success.

## Practical notes for this environment
- arxiv.org is NOT reachable via curl/WebFetch. Use the MCP tool
  `mcp__Exa__web_fetch_exa` (load it with ToolSearch "select:mcp__Exa__web_fetch_exa,mcp__Exa__web_search_exa")
  on URLs like https://arxiv.org/abs/<ID> (abstract) or https://arxiv.org/pdf/<ID> /
  https://arxiv.org/html/<ID> (full text, use maxCharacters 20000+ and search for equations).
  `mcp__Exa__web_search_exa` and `WebSearch` work for discovery.
- An arXiv ID counts as VERIFIED only if you fetched its abstract page in this session and the
  title/authors match what you claim. Otherwise mark it [UNVERIFIED-ID]. Never guess an ID.
- Equation numbers: give them when you actually read them in the fetched text; otherwise write
  "(eq. not checked)". Do not invent equation numbers.
- Python 3.13 with sympy, numpy, scipy, mpmath is installed. Save code under the directory
  you are told; run it; save outputs next to it.

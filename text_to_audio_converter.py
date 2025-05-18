import wave
import struct
import base64
import os

def text_to_wav(text_file, wav_file, bits_per_sample=16):
    """
    Convert a text file to a WAV file by encoding the text as binary data
    
    Args:
        text_file: Path to the input text file
        wav_file: Path to the output WAV file
        bits_per_sample: Audio quality (8, 16, or 32 bits)
    """
    # Read the text file
    with open(text_file, 'rb') as f:
        text_data = f.read()
    
    # Base64 encode the text data to ensure it can be properly represented
    encoded_data = base64.b64encode(text_data)
    
    # Add header to identify this as our encoded file
    header = b"ENCODED_TEXT:"
    data_with_header = header + encoded_data
    
    # Calculate number of audio samples needed
    num_samples = len(data_with_header)
    
    # Set audio parameters
    num_channels = 1  # Mono
    sample_rate = 44100  # CD quality
    
    # Create a new WAV file
    with wave.open(wav_file, 'w') as wav:
        wav.setnchannels(num_channels)
        wav.setsampwidth(bits_per_sample // 8)
        wav.setframerate(sample_rate)
        
        # Write each byte as an audio sample
        for byte in data_with_header:
            # Scale the value to fill the available bit range
            if bits_per_sample == 8:
                # 8-bit is unsigned bytes (0 to 255)
                wav.writeframes(struct.pack('B', byte))
            elif bits_per_sample == 16:
                # 16-bit uses signed integers (-32768 to 32767)
                # Scale from 0-255 to cover a reasonable range of the 16-bit space
                scaled_value = (byte - 128) * 256
                wav.writeframes(struct.pack('<h', scaled_value))
            elif bits_per_sample == 32:
                # 32-bit uses signed integers
                scaled_value = (byte - 128) * 16777216
                wav.writeframes(struct.pack('<i', scaled_value))
    
    print(f"Successfully converted {text_file} to {wav_file}")
    print(f"Original file size: {len(text_data)} bytes")
    print(f"Encoded WAV file size: {os.path.getsize(wav_file)} bytes")

def wav_to_text(wav_file, output_text_file):
    """
    Convert a WAV file back to the original text file
    
    Args:
        wav_file: Path to the input WAV file
        output_text_file: Path to the output text file
    """
    # Open the WAV file
    with wave.open(wav_file, 'r') as wav:
        num_channels = wav.getnchannels()
        sample_width = wav.getsampwidth()
        num_frames = wav.getnframes()
        
        # Read all frames
        raw_data = wav.readframes(num_frames)
        
        # Determine how to unpack the data based on sample width
        if sample_width == 1:
            # Process unsigned 8-bit audio
            values = []
            for i in range(0, len(raw_data), 1):
                values.append(struct.unpack('B', raw_data[i:i+1])[0])
            decoded_bytes = bytes(values)
        elif sample_width == 2:
            # Process signed 16-bit audio
            values = []
            for i in range(0, len(raw_data), 2):
                if i+2 <= len(raw_data):
                    value = struct.unpack('<h', raw_data[i:i+2])[0]
                    values.append(min(255, max(0, (value // 256) + 128)))
            decoded_bytes = bytes(values)
        elif sample_width == 4:
            # Process signed 32-bit audio
            values = []
            for i in range(0, len(raw_data), 4):
                if i+4 <= len(raw_data):
                    value = struct.unpack('<i', raw_data[i:i+4])[0]
                    values.append(min(255, max(0, (value // 16777216) + 128)))
            decoded_bytes = bytes(values)
    
    # Look for our header
    header = b"ENCODED_TEXT:"
    if header in decoded_bytes:
        # Extract the encoded data
        start_idx = decoded_bytes.find(header) + len(header)
        encoded_data = decoded_bytes[start_idx:]
        
        try:
            # Decode the Base64 encoded data
            original_data = base64.b64decode(encoded_data)
            
            # Write the original data to the output file
            with open(output_text_file, 'wb') as f:
                f.write(original_data)
            
            print(f"Successfully converted {wav_file} to {output_text_file}")
            print(f"Recovered file size: {len(original_data)} bytes")
        except Exception as e:
            print(f"Error decoding data: {e}")
    else:
        print("This WAV file doesn't contain encoded text data.")

# Example usage
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 3:
        print("Usage:")
        print("  To encode: python script.py encode input.txt output.wav")
        print("  To decode: python script.py decode input.wav output.txt")
        sys.exit(1)
        
    command = sys.argv[1].lower()
    
    if command == "encode":
        if len(sys.argv) < 4:
            print("Please provide input text file and output WAV file")
            sys.exit(1)
        text_to_wav(sys.argv[2], sys.argv[3])
    elif command == "decode":
        if len(sys.argv) < 4:
            print("Please provide input WAV file and output text file")
            sys.exit(1)
        wav_to_text(sys.argv[2], sys.argv[3])
    else:
        print("Unknown command. Use 'encode' or 'decode'")
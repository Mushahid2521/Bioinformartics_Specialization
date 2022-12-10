from bitstring import BitArray, BitStream

mapping = {'00': 'A', '01': 'T', '10': 'C', '11': 'G', 'A': '00', 'T': '01', 'C': '10', 'G': '11', '': ''}


def compress_genome(input_file):
    b = ''
    with open(input_file, 'r') as file:
        seq = file.readlines()[0]
        for c in seq:
            b += mapping[c]

    with open(input_file.split('.')[0], 'wb') as file:
        BitArray(bin=b).tofile(file)


def decompress_genome(input_byte, output_file):
    out_string = ""
    with open(input_byte, 'rb') as file:
        content = file.readlines()
        for chunk in content:
            bit_string = BitStream(chunk).bin
            for i in range(len(bit_string) + 1 // 2):
                out_string += mapping[bit_string[i * 2:i * 2 + 2]]

    with open(output_file, 'w') as file:
        file.write(out_string)


compress_genome('E_coli.txt')
decompress_genome('E_coli', "E_coli_output.txt")

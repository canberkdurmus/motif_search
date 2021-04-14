from dna import Dna
from motif_search import MotifSearch

if __name__ == '__main__':
    file_name = 'dna'
    dna = Dna(file_name, generate=False)
    dna_lines = dna.lines
    print('10-Mer: ', dna.ten_mer)
    print('Mutated 10-Mers: ', dna.mutated_list)

    for k in range(9, 12):
        MotifSearch('random', dna_lines, k, check_period=500).print_results()
        MotifSearch('gibbs', dna_lines, k, check_period=500).print_results()

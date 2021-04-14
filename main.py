from dna import Dna
from motif_search import MotifSearch

if __name__ == '__main__':
    file_name = 'dna'
    dna = Dna(file_name, generate=False).lines

    for k in range(9, 12):
        MotifSearch('random', dna, k, check_period=500).print_results()
        MotifSearch('gibbs', dna, k, check_period=500).print_results()

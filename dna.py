import random


class Dna:
    def __init__(self, file_name, generate=False):
        self.file_name = file_name
        self.lines = None
        self.ten_mer = None
        self.mutated_list = None
        if generate:
            self.generate_input_file()
        else:
            self.read_file()

    @staticmethod
    def get_random_nucleotide():
        index = random.randint(0, 3)
        nucleotides = ["a", "c", "g", "t"]
        return str(nucleotides[index])

    def read_file(self):
        lines = open(self.file_name + '.txt', "r").readlines()
        for line in lines:
            lines[lines.index(line)] = line.strip()

        info_list = open(self.file_name + '_info.txt', "r").readlines()
        self.ten_mer = info_list[0].strip()
        self.mutated_list = info_list[1].strip().split(' ')

        self.lines = lines
        return

    def mutate_k_mer(self, k_mer, num_of_mutation):
        mutated_positions = []
        for i in range(0, num_of_mutation):
            pos = random.randint(0, len(k_mer) - 1)
            while mutated_positions.__contains__(pos):
                pos = random.randint(0, len(k_mer) - 1)
            mutated_positions.append(pos)
            new_nuc = self.get_random_nucleotide()
            while k_mer[pos] == new_nuc:
                new_nuc = self.get_random_nucleotide()
            k_mer = k_mer[:pos] + new_nuc + k_mer[pos + 1:]
        return k_mer

    def generate_input_file(self):
        ten_mer = ''
        lines = []
        for i in range(0, 10):
            ten_mer += self.get_random_nucleotide()

        self.ten_mer = ten_mer
        self.mutated_list = []

        for i in range(0, 10):
            self.mutated_list.append(self.mutate_k_mer(ten_mer, 4))

        with open(str(self.file_name) + '.txt', 'w') as file:
            for i in range(0, 10):
                line = ""
                ten_mer_pos = random.randint(0, 490)
                for c in range(0, 500):
                    # Placing 10-mer
                    if c == ten_mer_pos:
                        line += str(self.mutated_list[i])
                    # Skip already filled positions
                    elif ten_mer_pos + 10 > c > ten_mer_pos:
                        continue
                    # Keep generating random nucleotides if no 10-mer placed
                    else:
                        line += self.get_random_nucleotide()
                lines.append(line)
                file.write(line + '\n')
        self.lines = lines
        with open(str(self.file_name) + '_info.txt', 'w') as file:
            file.write(ten_mer + '\n')
            file.write(str(self.mutated_list).replace('[', '').replace(']', '').replace('\'', '').replace(',', ''))
            file.write('\n')
        return

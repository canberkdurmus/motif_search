import random
import operator
import statistics


class MotifSearch:
    def __init__(self, search_type, dna, k, check_period=50):

        # Initial check for inputs
        if search_type != 'random' and search_type != 'gibbs':
            print('Invalid search type: ', search_type, ' should be \'gibbs\' or \'random\'')
            return

        # Initial object attitudes
        self.search_type = search_type
        self.dna = dna
        self.k = k
        self.check_period = check_period

        # Attitudes to use while making calculations
        self.current_profile = None
        self.temp_score_buffer = []

        # Final result attitudes
        self.final_score = None
        self.consensus_string = None
        self.best_motifs = None
        self.iterations = None

        # Primary jobs before running any motif search algorithm runs
        self.initial_motifs = self.select_initial_motifs()

        if self.search_type == 'random':
            self.randomized_motif_search()
        else:
            self.gibbs_sampler()

    def select_initial_motifs(self):
        motifs = []
        for dna_string in self.dna:
            pos = random.randint(0, len(dna_string) - self.k)
            motifs.append(dna_string[pos:pos + self.k])
        return motifs

    def set_new_profile(self, target_motifs):
        profile = [[], [], [], []]
        for i in range(0, len(target_motifs[0])):
            a = 0
            c = 0
            g = 0
            t = 0
            for motif in target_motifs:
                if motif[i] == "a":
                    a += 1
                elif motif[i] == "c":
                    c += 1
                elif motif[i] == "g":
                    g += 1
                elif motif[i] == "t":
                    t += 1
            if self.search_type == 'gibbs':
                a += 1
                c += 1
                g += 1
                t += 1
                total = len(target_motifs) + 4
            else:
                total = len(target_motifs)
            profile[0].append(a / total)
            profile[1].append(c / total)
            profile[2].append(g / total)
            profile[3].append(t / total)

        self.current_profile = profile
        return self.current_profile

    def calculate_prob_of_motif(self, motif):
        nucleotides_indices = {"a": 0, "c": 1, "g": 2, "t": 3}
        probability = 1
        for j in range(0, len(motif)):
            nuc = motif[j]
            try:
                probability *= self.current_profile[nucleotides_indices[nuc]][j]
            except ValueError:
                print(nuc, j)
        return probability

    def get_motif_from_string(self, dna_string, prob_dist=False):
        highest_prob = 0
        selected_motif = ""
        motif_prob_list = []
        for i in range(0, len(dna_string) - self.k):
            current_motif = dna_string[i:i + self.k]
            prob = self.calculate_prob_of_motif(current_motif)
            if prob_dist:
                motif_prob_list.append(prob)
            else:
                if prob > highest_prob:
                    highest_prob = prob
                    selected_motif = current_motif
        if prob_dist:
            indices = list(range(0, len(dna_string) - self.k))
            index = random.choices(indices, motif_prob_list)[
                0]  # Select a random number corresponding prob. distribution
            selected_motif = dna_string[index:index + self.k]
        return selected_motif

    def create_new_motifs(self):
        new_motifs = []
        for dna_string in self.dna:
            selected_motif = self.get_motif_from_string(dna_string)
            new_motifs.append(selected_motif)
        return new_motifs

    @staticmethod
    def calculate_score(target_motifs):
        score = 0
        for i in range(0, len(target_motifs[0])):
            nucleotide_counts = [0, 0, 0, 0]
            for motif in target_motifs:
                if motif[i] == "a":
                    nucleotide_counts[0] += 1
                elif motif[i] == "c":
                    nucleotide_counts[1] += 1
                elif motif[i] == "g":
                    nucleotide_counts[2] += 1
                elif motif[i] == "t":
                    nucleotide_counts[3] += 1
            score += len(target_motifs) - max(nucleotide_counts)
        return score

    def get_consensus_string(self):
        consensus = ""
        nucleotides = ["a", "c", "g", "t"]
        for i in range(0, len(self.best_motifs[0])):
            nuc_counts = [0, 0, 0, 0]
            for motif in self.best_motifs:
                if motif[i] == "a":
                    nuc_counts[0] += 1
                elif motif[i] == "c":
                    nuc_counts[1] += 1
                elif motif[i] == "g":
                    nuc_counts[2] += 1
                elif motif[i] == "t":
                    nuc_counts[3] += 1
            index, value = max(enumerate(nuc_counts), key=operator.itemgetter(1))
            consensus += str(nucleotides[index])
        return consensus

    def print_results(self):
        print('-------------------------------')
        print('Search Type: ', self.search_type, '  k=' + str(self.k))
        print('Final Score: ', self.final_score)
        print('Consensus String: ', self.consensus_string)
        print('Best Motifs: ', self.best_motifs)
        print('Number of Iterations: ', self.iterations)
        print('')

    def randomized_motif_search(self):
        motifs = self.initial_motifs
        iterations = 0
        while True:
            iterations += 1
            self.set_new_profile(motifs)
            motifs = self.create_new_motifs()
            score = self.calculate_score(motifs)
            self.temp_score_buffer.append(score)

            if len(self.temp_score_buffer) >= self.check_period:
                mean = statistics.mean(self.temp_score_buffer)
                # print('== Random Check ==     mean: ', mean, 'score: ', score)
                if mean <= score:
                    # Terminate random search
                    self.best_motifs = motifs
                    self.iterations = iterations
                    self.final_score = score
                    self.consensus_string = self.get_consensus_string()
                    return
                else:
                    self.temp_score_buffer.clear()

    def gibbs_sampler(self):
        motifs = self.initial_motifs
        iterations = 0
        while True:
            iterations += 1

            # Select a random motif from the list
            i = random.randint(0, len(self.dna) - 1)
            copy_motifs = motifs.copy()
            copy_motifs.pop(i)

            # Profile motifs except selected one
            self.set_new_profile(copy_motifs)

            # Add new motif to empty place
            motifs[i] = self.get_motif_from_string(self.dna[i], prob_dist=True)
            score = self.calculate_score(motifs)
            self.temp_score_buffer.append(score)

            if len(self.temp_score_buffer) >= self.check_period:
                mean = statistics.mean(self.temp_score_buffer)
                # print('== Gibbs Check ==     mean: ', mean, 'score: ', score)
                if mean <= score:
                    # Terminate random search
                    self.best_motifs = motifs
                    self.iterations = iterations
                    self.final_score = score
                    self.consensus_string = self.get_consensus_string()
                    return
                else:
                    self.temp_score_buffer.clear()

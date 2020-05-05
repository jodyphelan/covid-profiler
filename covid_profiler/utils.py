from collections import Counter

def get_N_content(seq):
    cnt = Counter(seq)
    return cnt["N"]/len(seq)

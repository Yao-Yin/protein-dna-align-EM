import random as rd
valid_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\pg_500_hg19_presuf10_valid.txt"
training_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\pg_500_hg19_presuf10_train.txt"
test_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\pg_500_hg19_presuf10_test.txt"
ori_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\pg_500_hg19_presuf10.txt"

rd.seed(10)
removeRepeats = set()

def write_data(idxs, seq, filename):
    with open(filename, "w") as f:
        #get 4 seqs
        for i in idxs:
            curr_dna = seq[4*i+3].rstrip()
            curr_pro = seq[4*i+1].rstrip()
            if (curr_pro, curr_dna) not in removeRepeats:
                f.write(seq[4*i])
                f.write(seq[4*i+1])
                f.write(seq[4*i+2])
                f.write(seq[4*i+3].rstrip() + "\n")
                removeRepeats.add((curr_pro, curr_dna))



with open(ori_file, "r") as f:
    all_seqs = f.readlines()
    all_seqs[-1] += "\n"
    tot = len(all_seqs) // 4
    print(tot)
    idx = [i for i in range(tot)]
    rd.shuffle(idx)
    train_num = int(tot*0.6)
    valid_num = int(tot*0.2)
    test_num = tot-train_num-valid_num
    train_idx = [idx[i] for i in range(train_num)]
    valid_idx = [idx[i] for i in range(train_num, train_num + valid_num)]
    test_idx = [idx[i] for i in range(train_num+valid_num, tot)]
    print(len(train_idx))
    print(len(valid_idx))
    print(len(test_idx))
    write_data(train_idx, all_seqs, training_file)
    write_data(valid_idx, all_seqs, valid_file)
    write_data(test_idx, all_seqs, test_file)
    train_idx.extend(valid_idx)
    train_idx.extend(test_idx)
    train_idx.sort()
    for i in range(len(train_idx)):
        if train_idx[i] != i:
            print(i)
    

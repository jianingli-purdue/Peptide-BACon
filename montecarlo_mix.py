from grid import grid
import dev
from box import box
import sys
import os
import random
from math import exp
from scipy.constants import Boltzmann

# MONTECARLO ALGORITHM FOR CO-ASSEMBLY
##------------------------------------------------------------------------------------------------------------
# write sequences to steps file
def write_seqs(current_seqs, flag):
    steps_file = open("montecarlo/" + name + "_steps.txt", "a")
    i = 0
    for seq in current_seqs:
        i += 1
        steps_file.write(str(seq))
        if i != len(current_ses): # check to make sure this works
            steps_file.write(",")
    if flag:
        steps_file.write(":" + str(current_pep_AP) + ":y\n") # accepted
    elif not flag:
        steps_file.write(":" + str(current_pep_AP) + ":n\n") # rejected
    steps_file.close()

##------------------------------------------------------------------------------------------------------------
# function to generate name from list of peptides
def gen_name(peps):
    name = ""
    for pep in peps: # iterate through all peptides
        for aa in pep:  # iterate through peptide's amino acids
            name += str(aa)
        name += "_"
    return name

##------------------------------------------------------------------------------------------------------------
# function to build all possible new peptides from potential group
def groupset(peps, group, position=None, pep=None):
    set = []
    if position is None:
        for res in group:
            tmp = []
            # add peptide we are not modifying to tmp
            for i in range(2):
                if i != position:
                    tmp.append(peps[i])
            # loop and create modified pairs
            for to_mutate in range(len(peps[pep])):
                new_pep = ""
                for i in range(len(peps[pep])): # big big ohs only
                    if i == to_mutate:
                        new_pep += res
                    else:
                        new_pep += peps[pep][i]
                tmp.append(new_pep)
                set.append(tmp)
                set.append(make_reverse(tmp))
    else:
        # initialize temporary variables
        to_mutate = position
        # loop through potential mutations
        for res in group:
            tmp = []
            # add peptide we are not modifying to tmp
            for i in range(2):
                if i != pep:
                    tmp.append(peps[i])
            # add modified peptide to tmp
            new_pep = ""
            for i in range(len(peps[pep])): # big big ohs only
                if i == to_mutate:
                    new_pep += res
                else:
                    new_pep += peps[pep][i]
            tmp.append(new_pep)
            set.append(tmp)
            set.append(make_reverse(tmp))
    return set

##------------------------------------------------------------------------------------------------------------
# function to make complementary list
def make_reverse(list):
    return list[::-1]

##------------------------------------------------------------------------------------------------------------
# function to find group of amino acid
def find_group(aa):
    for group in aa_groups:
        for res in group:
            if aa == res:

                return group

##------------------------------------------------------------------------------------------------------------
# function to check for duplicates in list of previously tried peptides
def check_duplicates(new_peps, prev_tried):
    # create reversed list, as position of peptides is irrelevant
    tmp = [new_peps[-1], new_peps[0]]
    # print("checking duplicates")
    # print(new_peps)
    if (new_peps in prev_tried):
        # print("true")
        return True
    elif (tmp in prev_tried):
        # print("true")
        return True
    else:
        # print("false")
        return False

##------------------------------------------------------------------------------------------------------------
# function to flatten list of lists
def flatten(lists):
    if lists == None:
        return None
    to_return = []
    for list in lists:
        for item in list:
            to_return.append(item)
    return to_return

##------------------------------------------------------------------------------------------------------------
# function to copy list
def copy_list(list):
    to_return = []
    for item in list:
        to_return.append(item)
    return to_return

##------------------------------------------------------------------------------------------------------------
# generate new ste
def gen_step(old_peps, old_groups, accepted, position=None, pep=None, flag=False):
    # if we did not accept previous mutation
    if not accepted:
        # flag indicates that we should keep generating mutations from the same group
        # if all mutations at this pep/position have been made, set flag false

        all_mutations = groupset(old_peps, find_group(old_peps[pep][position]), position=position, pep=pep)
        tmp = []
        for mutation in all_mutations:
            tmp.append(check_duplicates(mutation, prev_tried))
        # if all mutations from group have been attempted, set flag False
        if all(tmp):
            flag = False

        # if all previous groups have been tested then start with a random group
        if len(old_groups) >= 6 and flag is False:
            new_group = aa_groups[random.randint(0, 5)]
            # test whether all permutations of new group have already been tried
            change_group = True
            while change_group is True:
                new_group = aa_groups[random.randint(0,5)]
                all_mutations = groupset(old_peps, new_group, pep=pep)
                # check to see if we've made all mutations from new group
                tmp = []
                for mutation in all_mutations:
                    tmp.append(check_duplicates(mutation, prev_tried))
                if not all(tmp):
                    change_group = False
        # otherwise pick an untested group
        else:
            if flag:
                return gen_step(old_peps, None, True, position=position, pep=pep, flag=flag)
            else:
                new_group = old_groups[0]
                while new_group in old_groups:
                    new_group = aa_groups[random.randint(0, 5)]

        # generate new pep
        new_peps = copy_list(old_peps)
        seq_length = len(new_peps[pep])

        while new_peps == old_peps or check_duplicates(new_peps, prev_tried):
            posn_to_mutate = random.randint(0, seq_length)
            pep_to_mutate = random.randint(0, 1)
            res = new_group[random.randint(0, len(new_group)-1)]
            new_pep = ""
            for i in range(len(old_peps[pep_to_mutate])):
                if i == posn_to_mutate:
                    new_pep += res
                else:
                    new_pep += old_peps[pep_to_mutate][i]
            #print(new_pep)
            new_peps[pep_to_mutate] = new_pep

    # if we did accept previous mutation
    elif accepted:
        # keep group
        group = find_group(old_peps[pep][position])
        # marker to check if we have tried all substitutions from one group
        attempted_aas = []
        # generate new pep
        new_peps = copy_list(old_peps)

        while check_duplicates(new_peps, prev_tried) and (set(attempted_aas) != set(group)):
            posn_to_mutate = position
            pep_to_mutate = pep
            res = group[random.randint(0, len(group))-1]
            while res in attempted_aas: # make sure we're not trying the same res twice
                res = group[random.randint(0, len(group))-1]
            attempted_aas.append(res)
            new_pep = ""
            for i in range(len(old_peps[pep_to_mutate])):
                if i == posn_to_mutate:
                    new_pep += res
                else:
                    new_pep += old_peps[pep_to_mutate][i]
            new_peps[pep_to_mutate] = new_pep

        # if we have tried all amino acids of current group go to fresh gen
        if check_duplicates(new_peps, prev_tried) and set(attempted_aas) == set(group):
            new_peps, posn_to_mutate, pep_to_mutate, flag = gen_step(old_peps, group, False, position=position, pep=pep_to_mutate, flag=flag)

    return new_peps, posn_to_mutate, pep_to_mutate, flag

##------------------------------------------------------------------------------------------------------------
## compare potential step to prev current
def compare(new_AP, old_AP, step):
    delta_AP = new_AP - old_AP
    # metropolis function - initally accepts about 1/3 of positive delta_APs
    total_steps = 150
    if step > 40:
        addn = 3.5
    else:
        addn = 3.5/total_steps * step
    bolt_val = 1/(1+exp(delta_AP - 0.5 + addn))
    # check if downhill move
    if (delta_AP) <= 0:
        return True
    # metropolis function used to accept uphill moves
    elif bolt_val >= random.random():
        return True
    else:
        return False

##------------------------------------------------------------------------------------------------------------
## function to find filepath
def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

##------------------------------------------------------------------------------------------------------------
# global lists of amino acids
amino_acids = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
positive_charge_aa = ["R", "H", "K"]
negative_charge_aa = ["D", "E"]
polar_uncharged_aa = ["S", "T", "N", "Q", "C"]
special_aa = ["G", "P"]
hydrophobic_aa = ["A", "V", "I", "L", "M"]
hydrophobic_ring_aa =  ["F", "Y", "W"]
aa_groups = [positive_charge_aa, negative_charge_aa, polar_uncharged_aa, special_aa, hydrophobic_aa, hydrophobic_ring_aa]

# initialize attempt and path
attempt = 0
path = "/home/marlo/HDD/Maestro/jobs/montecarlo/"

# create lists of tried and accepted sequences
prev_tried = []
steps = []
flag = False

# get trial_name of run for record and make directory if necessary
#trial_name = str(input("Trial name: "))
trial_name = input("Trialname: ")

try:
    os.mkdir(path + trial_name + '/')
    shutil.move(find("mr_md.msj", "home/"), path + trial_name + "/mr_md.msj")
    shutil.move(find("mr_minimization.msj", "home/"), path + trial_name + "/mr_minimization.msj")
except:
    pass

# check to see if file with previous steps exists
exists = os.path.isfile("montecarlo/" + trial_name + "_steps.txt")

# if previous steps file exists, load in all previous attempts
if exists:
    # specificy sequence to start with in command prompt (MUST HAVE ALREADY BEEN SIMULATED)
    current_seq1, current_seq2 = sys.argv[1], sys.argv[2]
    current_seqs = [current_seq1, current_seq2]
    step = -1
    # load list of attempts
    steps_file = open("montecarlo/" + trial_name + "_steps.txt", "r")
    for line in steps_file:
        temp = line.split(":")
        prev_tried.append(list(temp[0])) # this may not work
        # update attempt and step number
        if temp[2].strip() == "y":
            step += 1
            steps.append([current_seq1, current_seq2])
            if temp[0].strip() == (current_seq1 + "," + current_seq2):
                init_step = step
                init_attempt = 2
            attempt = 0
        elif temp[2].strip() == "n":
            attempt += 1
    steps_file.close()

    # analyze specified current sequence
    current_peps = box(name=gen_name(current_seqs) + "step" + str(init_step) + "_attempt" + str(init_attempt), spacing=15, height=8, total=512, rand=False, aaCG=True, trialName=trial_name, path=path)
    current_AP = current_peps.analyze()
    #current_AP = random.random()

    # correct attempts and steps
    if attempt != 0:
        attempt += 1
    else:
        step += 1
else:
    # generate random sequence
    step = 0
    #seq_length1 = int(input("First sequence length: "))
    seq_length1 = 2
    current_seq1 = ""
    for i in range(seq_length1):
        current_seq1 += amino_acids[random.randint(0, 19)]
    #seq_length2 = int(input("First sequence length: "))
    seq_length2 = 2
    current_seq2 = ""
    for i in range(seq_length2):
        current_seq2 += amino_acids[random.randint(0, 19)]
    current_seqs = [current_seq1, current_seq2]

    print(current_seqs)
    # build and run first box
    current_peps = box(name=gen_name(current_seqs) + "step" + str(step) + "_attempt" + str(attempt), spacing=15, height=8, total=512, rand=False, aaCG=True, trialName=trial_name , path=path)
    current_peps.addStruct(current_seqs, list=True)
    current_peps.build()
    current_peps.sysBuild()
    current_peps.run()
    current_AP = current_peps.analyze()
    #current_AP = random.random()

    # append to list of steps and previously tried
    steps.append(current_seqs)
    prev_tried.append(current_seqs)
    steps_file = open("montecarlo/" + trial_name + "_steps.txt", "a")
    steps_file.write(str(current_seqs[0]) + "," + str(current_seqs[1]) + ":" + str(current_AP) + ":y" + "\n")
    steps_file.close()

print("Starting Monte Carlo method!")
print("Seq: " + current_seq1)
print("Seq: " + current_seq2)
print("AP:" + str(current_AP))

# generate first new peptide
new_job = current_peps
new_seqs = [current_seq1, current_seq2]
pep_to_mutate = random.randint(0,1)
posn_to_mutate = random.randint(0,1)
new_AP = current_AP
step += 1

# get amino acid group of first new peptide
tried_groups = [find_group(new_seqs[pep_to_mutate][posn_to_mutate])]

# take steps
for i in range(150):
    # append previous sequence to prev_tried
    prev_tried.append(new_seqs)
    print(steps)

    # run sysbuilder in case we reject current trial
    rejected_seqs, rejected_posn, rejected_pep, flag = gen_step(steps[-1], tried_groups, False, position=posn_to_mutate, pep=pep_to_mutate, flag=flag)
    rejected_group = find_group(rejected_seqs[rejected_pep][rejected_posn])
    rej_job = box(name=gen_name(rejected_seqs) + "step" + str(step) + "_attempt" + str(attempt+1), spacing=15, height=8, total=512, rand=True, aaCG=True, concentrations=[0.5, 0.5], trialName=trial_name, path=path)
    print(rejected_seqs)
    for seq in rejected_seqs:
        rej_job.addStruct(seq)
    rej_job.build()
    rej_job.sysBuild()
    print("potential rej: " + str(rejected_seqs))

    # run sysbuilder in case we accept current trial
    accepted_seqs, accepted_posn, accepted_pep, flag = gen_step(new_seqs, None, True, position=posn_to_mutate, pep=pep_to_mutate, flag=flag)
    accepted_group = find_group(accepted_seqs[accepted_pep][accepted_posn])
    acc_job = box(name=gen_name(accepted_seqs) + "step" + str(step+1) + "_attempt0", spacing=15, height=8, total=512, rand=True, aaCG=True, concentrations=[0.5, 0.5], trialName=trial_name, path=path)
    for seq in accepted_seqs:
        acc_job.addStruct(seq)
    acc_job.build()
    acc_job.sysBuild()
    print("potential acc: " + str(accepted_seqs))


    # analyze latest run
    new_AP = new_job.analyze()
    #new_AP = random.random()

    print("Seq: " + gen_name(new_seqs))
    print("AP: " + str(new_AP))
    # compare new to current
    if compare(new_AP, current_AP, step):
        # take step
        steps.append(new_seqs)
        current_seqs = new_seqs
        current_job = new_job
        current_AP = new_AP
        step += 1
        attempt = 0

        # append to list of steps
        steps_file = open("montecarlo/" + trial_name + "_steps.txt", "a")
        steps_file.write(str(current_seqs[0]) + "," + str(current_seqs[1]) + ":" + str(current_AP) + ":y" + "\n")
        steps_file.close()
        print("accepted")

        # start new job
        new_job = acc_job
        new_seqs = accepted_seqs
        posn_to_mutate = accepted_posn
        pep_to_mutate = accepted_pep
        tried_groups = [find_group(new_seqs[pep_to_mutate][posn_to_mutate])]
        flag = True

        print("Trying: " + str(new_seqs) + "\n")
        new_job.run()
        rej_job.deleteJob()
    else:
        # append list of steps
        steps_file = open("montecarlo/" + trial_name + "_steps.txt", "a")
        steps_file.write(str(new_seqs[0]) + "," + str(new_seqs[1]) + ":" + str(new_AP) + ":n" + "\n")
        steps_file.close()
        attempt += 1
        tried_groups.append(find_group(new_seqs[pep_to_mutate][posn_to_mutate]))
        print("rejected")

        # start new job
        new_job = rej_job
        new_seqs = rejected_seqs
        posn_to_mutate = rejected_posn
        pep_to_mutate = rejected_pep
        print("Trying: " + str(new_seqs) + '\n')
        new_job.run()
        acc_job.deleteJob()

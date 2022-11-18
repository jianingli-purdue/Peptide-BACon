import os
import dev
import sys
import shutil
import random
from box import box
from math import exp
from grid import grid
# from scipy.constants import Boltzmann
## MONTECARLO ALGORITHM FOR SINGLE PEPTIDES
##------------------------------------------------------------------------------------------------------------
# function to build all possible new peptides from potential group
def groupset(pep, group, position=None):
    set = []
    if position is None:
        for new_res in group:
            for to_mutate in range(len(pep)):
                new_pep = ""
                for i in range(len(pep)): # big big ohs only
                    if i == to_mutate:
                        new_pep += new_res
                    else:
                        new_pep += pep[i]
                set.append(new_pep)
    else:
        to_mutate = position
        for new_res in group:
            new_pep = ""
            for i in range(len(pep)): # big big ohs only
                if i == to_mutate:
                    new_pep += new_res
                else:
                    new_pep += pep[i]
            set.append(new_pep)
    return set

##------------------------------------------------------------------------------------------------------------
# function to find group of amino acid
def find_group(aa):
    for group in aa_groups:
        for res in group:
            if aa == res:
                return group

##------------------------------------------------------------------------------------------------------------
# generate new step
def gen_step(old_pep, old_groups, accepted, position=None, flag=False):
    # if we did not accept previous mutation
    if not accepted:
        # flag indicates that we should keep generating mutations from the same group
        # if all possible mutations from a group is a subset of previously tried mutations, set flag False
        if set(groupset(old_pep, find_group(old_pep[position]), position=position)).issubset(set(prev_tried)):
            flag = False

        # if all groups have been tried and flag is false, pick a group randomly
        if len(old_groups) >= 6 and flag is False:
            new_group = aa_groups[random.randint(0, 5)]
            # test whether all permutations of new group have already been tried
            while set(groupset(old_pep, new_group)).issubset(set(prev_tried)):
                new_group = aa_groups[random.randint(0, 5)]

        # otherwise either pick from same group if flag is True or pick random untried group
        else:
            if flag:
                return gen_step(old_pep, None, True, position=position, flag=flag)
            else:
                new_group = old_groups[0]
                while new_group in old_groups:
                    new_group = aa_groups[random.randint(0, 5)]

        # generate new pep
        new_pep = old_pep
        seq_length = len(new_pep)

        # while new pep is same as old or in prev tried, generate new pep
        # should be able to change this to just "new_pep in prev_tried" (?)
        while new_pep == old_pep or new_pep in prev_tried:
            to_mutate = random.randint(0, seq_length)
            new_res = new_group[random.randint(0, len(new_group)-1)]
            new_pep = ""
            for i in range(len(old_pep)):
                if i == to_mutate:
                    new_pep += new_res
                else:
                    new_pep += old_pep[i]

    # if we did accept previous mutation
    elif accepted:
        # keep group
        group = find_group(old_pep[position])

        # marker to check if we have tried all substitutions from one group
        attempted_aas = []

        # generate new pep
        new_pep = old_pep
        while (new_pep in prev_tried) and (set(attempted_aas) != set(group)):
            to_mutate = position
            new_res = group[random.randint(0, len(group))-1]
            while new_res in attempted_aas: # make sure we're not trying the same res twice
                new_res = group[random.randint(0, len(group))-1]
            attempted_aas.append(new_res)
            new_pep = ""
            for i in range(len(old_pep)):
                if i == to_mutate:
                    new_pep += new_res
                else:
                    new_pep += old_pep[i]

        # if we have tried all amino acids of current group go to fresh gen
        if (new_pep in prev_tried) and set(attempted_aas) == set(group):
            new_pep, to_mutate, flag = gen_step(old_pep, group, False, position=position, flag=flag)

    return new_pep, to_mutate, flag

##------------------------------------------------------------------------------------------------------------
## compare potential step to prev current
def compare(new_sasa, old_sasa, step):
    delta_sasa = new_sasa - old_sasa
    # metropolis function - initally accepts about 1/3 of positive delta_sasas
    total_steps = 150
    if step > 40:
        addn = 3.5
    else:
        addn = 3.5/total_steps * step
    bolt_val = 1/(1+exp(delta_sasa - 0.5 + addn))
    # check if downhill move
    if (delta_sasa) <= 0:
        return True
    # metropolis function used to accept uphill moves
    elif bolt_val >= random.random():msys.Save(msys_model, traj_location_ext + new_name + "pdb", structure_only=True)

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
# global lists of amino acidsG
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
trial_name = str(input("Trial name: "))
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
    current_seq = sys.argv[1]
    step = 11
    # load list of attempts
    steps_file = open("montecarlo/" + trial_name + "_steps.txt", "r")
    for line in steps_file:
        temp = line.split(":")
        prev_tried.append(temp[0])
        # update attempt and step number
        if temp[2].strip() == "y":
            step += 1
            steps.append(temp[0])
            if temp[0].strip() == current_seq:
                print("found starting seq")
                init_step = step
                init_attempt = attempt
            attempt = 0
        elif temp[2].strip() == "n":
            attempt += 1
    steps_file.close()
    # analyze specified current sequence
    current_pep = box(name=current_seq + "_step" + str(init_step) + "_attempt" + str(init_attempt), spacing=15, height=4, total=64, rand=False, aaCG=False, trialName=trial_name, path=path)
    current_pep_sasa = current_pep.analyze()


    # correct attempts and steps
    if attempt != 0:
        attempt += 1
    else:
        step += 1
else:
    # generate random sequence
    step = 0
    seq_length = int(input("Sequence length: "))
    current_seq = ""
    for i in range(seq_length):
        current_seq += amino_acids[random.randint(0, 19)]

    # build and run first box
    current_pep = box(name=current_seq + "_step" + str(step) + "_attempt" + str(attempt), spacing=15, height=4, total=64, rand=False, aaCG=False, trialName=trial_name , path=path)
    current_pep.addStruct(current_seq)
    current_pep.build()
    current_pep.sysBuild()
    current_pep.run()
    current_pep_sasa = current_pep.analyze()

    # append to list of steps and previously tried
    # steps_file = open("montecarlo/" + trial_name + "_steps.txt", "a")
    # steps_file.write(str(current_seq) + ":" + str(current_pep_sasa) + ":start" + "\n")
    # steps_file.close()
    steps.append(current_seq)
    prev_tried.append(current_seq)

print("Starting Monte Carlo method!")
print("Seq: " + current_seq)
print("SASA:" + str(current_pep_sasa))

# generate first new peptide
new_pep = current_pep
new_seq = current_seq
new_posn = random.randint(0,1)
new_pep_sasa = current_pep_sasa

step += 1


# get amino acid group of first new peptide
tried_groups = [find_group(new_seq[new_posn])]


# take steps
for i in range(200):
    # append previous sequence to prev_tried
    prev_tried.append(new_seq)

    # run sysbuilder in case we reject current trial
    new_rejected_seq, new_rejected_posn, flag = gen_step(steps[-1], tried_groups, False, position=new_posn, flag=flag)
    new_rejected_group = find_group(new_rejected_seq[new_rejected_posn])
    rej_pep = box(name=new_rejected_seq + "_step" + str(step) + "_attempt" + str(attempt+1), spacing=15, height=4, total=64, rand=False, aaCG=False, trialName=trial_name, path=path)
    rej_pep.addStruct(new_rejected_seq)
    rej_pep.build()
    rej_pep.sysBuild()
    print("potential rej: " + new_rejected_seq)


    # run sysbuilder in case we accept current trial
    new_accepted_seq, new_accepted_posn, flag = gen_step(new_seq, None, True, position=new_posn, flag=flag)
    new_accepted_group = find_group(new_accepted_seq[new_accepted_posn])
    acc_pep = box(name=new_accepted_seq + "_step" + str(step+1) + "_attempt0", spacing=15, height=4, total=64, rand=False, aaCG=False, trialName=trial_name, path=path)
    acc_pep.addStruct(new_accepted_seq)
    acc_pep.build()
    acc_pep.sysBuild()
    print("potential acc: " + new_accepted_seq)

    # analyze latest run
    new_pep_sasa = new_pep.analyze()
    print("Seq: " + new_seq)
    print("Score: " + str(1/new_pep_sasa))
    # compare new to current
    if compare(new_pep_sasa, current_pep_sasa, step):
        # take step
        steps.append(new_seq)
        current_seq = new_seq
        current_pep = new_pep
        current_pep_sasa = new_pep_sasa
        step += 1
        attempt = 0
        # append to list of steps
        if i > 1:
            steps_file = open("montecarlo/" + trial_name + "_steps.txt", "a")
            steps_file.write(str(current_seq) + ":" + str(current_pep_sasa) + ":y" + "\n")
            steps_file.close()
            print("accepted")
        # start new job
        new_pep = acc_pep
        new_seq = new_accepted_seq
        new_posn = new_accepted_posn
        tried_groups = [find_group(new_seq[new_posn])]
        flag = True

        print("Trying: " + new_seq + "\n")
        new_pep.run()
        # rej_pep.deleteJob()
    else:
        # append list of steps
        steps_file = open("montecarlo/" + trial_name + "_steps.txt", "a")
        steps_file.write(str(new_seq) + ":" + str(new_pep_sasa) + ":n" + "\n")
        steps_file.close()
        attempt += 1
        tried_groups.append(find_group(new_seq[new_posn]))
        print("rejected")
        # start new job
        new_pep = rej_pep
        new_seq = new_rejected_seq
        new_posn = new_rejected_posn
        print("Trying: " + new_seq)
        new_pep.run()
        # acc_pep.deleteJob()

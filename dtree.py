# Decision Tree Model for Data Set

# imports here
from numpy.core.fromnumeric import sort
from numpy.lib.arraysetops import unique
import pandas as pd
import numpy as np
from collections import Counter
import math
import os
import shutil

######################  functions to model the task  ######################


# function to calculate entropy of a particular dataset
def GetEntropy(df):
    N = len(df)
    hashMap = Counter(df["Output"])

    entropy = 0
    for item, val in hashMap.items():
        entropy -= (val / N) * math.log2((val / N))

    return entropy


def summary(case, treeTraversal, significanceList, EntropyList):
    info = ""
    file = open(path1, 'a')
    for i in range(len(EntropyList)):
        if i == 0:
            info += "Raw Data : [" + str(EntropyList[i]) + "]"
        else:
            info += " -> " + treeTraversal[i - 1] + \
                " : [" + str(EntropyList[i]) + "]"
    info += "\n\n"
    file.write(info)
    file.close()

    info = ""
    file = open(path2, 'a')
    for i in range(len(significanceList)):
        if i != 0:
            info += " -> "
        info += treeTraversal[i] + \
            " : [" + str(EntropyList[i]) + "]"
    info += "\n\n"
    file.write(info)
    file.close()

    for i in treeTraversal:
        print(i, end=" ")
        if i != treeTraversal[-1]:
            print("->", end=" ")
    if case == 0:
        print(" ...Splitted till the End!")
    elif case == 1:
        print(" [S <= epsilon (S)]")
    else:
        print(" [Len <= epsilon (Len)]")


# recursive function for tree modelling
def dTree(df, usedDecisions, treeTraversal, significanceList, EntropyList, path):

    df.to_excel(str(path + "/dataset.xlsx"))
    entropyS = GetEntropy(df)
    N = len(df)

    if len(treeTraversal) == 0:
        EntropyList.append(entropyInitial)

    file = open(str(path + "/info.txt"), 'w')
    info = _str1(N, entropyS)
    info += Probabilities_str(df)

    if usedDecisions.count(0) == 0:
        info += baseCond_str(0)
        file.write(info)
        file.close()

        summary(0, treeTraversal, significanceList, EntropyList)
        return

    if entropyS <= epsilonEntropy:
        info += baseCond_str(1)
        file.write(info)
        file.close()

        summary(1, treeTraversal, significanceList, EntropyList)
        return

    if N <= epsilonLen:
        info += baseCond_str(2)
        file.write(info)
        file.close()

        summary(2, treeTraversal, significanceList, EntropyList)
        return

    info += _str2()

    infoGainList = [0] * numAttributes

    for i in range(numAttributes):
        # A -> already splitted in
        if usedDecisions[i] == 1:
            info += alreadySplit_str(Enums[i])
            continue

        info += f"{Enums[i]} : "
        entropyCurrent = [0] * numVariants
        for j in range(numVariants):
            df_splitted = df.loc[df[Enums[i]] == variants[j]]
            entropyCurrent[j] = [GetEntropy(df_splitted), len(df_splitted)]

        # calculating information gain for this Attribute A
        ig = entropyS

        # split after : A1 A2 A3 A4 A5
        for entropy, length in entropyCurrent:
            info += f"{entropy}\t"
            if length > 0:
                ig -= (length / N) * entropy

        info += "\n"

        # information gain calculated for attribute Suppose A
        infoGainList[i] = ig

    info += _str3()

    max_ig = -1 * math.inf
    maxIg_attr_idx = -1

    for i in range(numAttributes):
        info += f"{Enums[i]} : {infoGainList[i]}\n"
        if max_ig < infoGainList[i] and usedDecisions[i] == 0:
            max_ig = infoGainList[i]
            maxIg_attr_idx = i

    # Ig of E being MAX :: E is to be splitted next
    usedDecisions[maxIg_attr_idx] = 1
    toSplitEnum = Enums[maxIg_attr_idx]

    info += determiningAttr_str(toSplitEnum, max_ig)

    # E1 E2
    for i in range(numVariants):
        df_splitted = df.loc[df[toSplitEnum] == variants[i]]

        info += lenAfterSplit_str(toSplitEnum, i, df_splitted)
        info += Significance_str(max_ig, len(df_splitted), numDataSet)

        # modify path
        path += f"/{toSplitEnum}{i + 1}"

        treeTraversal.append(str(toSplitEnum + str(variants[i])))
        significanceList.append(max_ig * (len(df_splitted) / numDataSet))
        EntropyList.append(GetEntropy(df_splitted))

        os.mkdir(path)
        dTree(df_splitted, usedDecisions, treeTraversal,
              significanceList, EntropyList, path)

        path = path.rsplit('/', 1)[0]
        treeTraversal.pop()
        significanceList.pop()
        EntropyList.pop()

    # write to the info.txt file
    file.write(info)
    file.close()

    # backtrack
    usedDecisions[maxIg_attr_idx] = 0


# Main function for Testing
def main():
    filePath = path = os.path.join(os.getcwd(), 'Inp1.xlsx')
    df = pd.read_excel(filePath)

    # To be taken as Input after words -> 0 && 1
    global nclassEntropy
    nclassEntropy = 2

    c = Counter(df['Output'])
    print(c.most_common(10))

    global possibleOutcomes
    possibleOutcomes = list(Counter(df['Output']).keys())
    possibleOutcomes.sort()
    nclassEntropy = len(Counter(df["Output"]))

    # Corresponds to A B C D E
    global numAttributes
    numAttributes = 2
    # numAttributes = int(input("Enter the Number of Attributes : "))

    # Corresponds to 5 in example
    global numVariants
    numVariants = 2
    global variants 
    variants = []
    for col in df:
        if col.lower()[:5] != "input" and col.lower()[:6] != 'output':
            variants += list(df[col].unique())
    variants = unique(variants)
    numVariants = len(variants)
    print("Variants : " + str(variants) , end="\n\n")
    # numVariants = int(input("Enter the Number of Variances : "))

    # Arbitarary
    global numDataSet
    numDataSet = len(df)

    global entropyInitial
    entropyInitial = GetEntropy(df)

    # Names A,B,C,D,E,.....
    global Enums
    Enums = list(df.columns[2:])
    print(Enums)

    # Epsilon Value for Length of DataSet to Split and Recurse
    global epsilonLen
    epsilonLen = 0
    # epsilonLen = int(input("Enter the Epsilon Value epsilon for Length of DataSet : "))

    # Epsilon Value for Entropy of DataSet to Split and Recurse
    global epsilonEntropy
    epsilonEntropy = 0
    # epsilonEntropy = int(input("Enter the Epsilon Value epsilon for Entropy of DataSet : "))

    # make the Results Folder for Tree Structure
    path = os.path.join(os.getcwd(), 'Results')

    global path1
    path1 = os.path.join(os.getcwd(), 'Summary1')
    fp = open(path1, "w")
    fp.write("Decision Modelling Summary Based on Entropy Changes.\n\n")
    fp.close()

    global path2
    path2 = os.path.join(os.getcwd(), 'Summary2')
    fp = open(path2, "w")
    fp.write("Decision Modelling Summary Based on Significance of Decision.\n\n")
    fp.close()

    if os.path.exists(path) and os.path.isdir(path):
        shutil.rmtree(path)

    os.mkdir(path)

    usedDecisions = [0] * numAttributes
    dTree(df, usedDecisions, [], [], [], path)


#############  String functions Useful in Printing Info.txt  #############


def _str1(N, entropyS):
    return f"PART - I\nInformation About the DataSet:\n\n- Number of Points in DataSet = {N}\n- Number of Attributes = {numAttributes}\n- Number of Variants = {numVariants}\n\n- Entropy(Data) = {entropyS}\n\nProbability of Occurence of Each Outcome in the DataSet is:\n"


def _str2():
    return f"\n\nPART - II\nInformation about Processing DataSet:\n\nEntropy Matrix: (attributes x variants)\n\n"


def _str3():
    return "\nInformation Gain for Each Attribute after Processing...\n"


def alreadySplit_str(enum):
    return f"{enum} : already splitted.\n"


def determiningAttr_str(toSplitEnum, max_ig):
    return f"\nDetermining Attribute for this DataSet is: {toSplitEnum}\nMaximum Information Gain Value of {max_ig}\n\nSplitting Data on the Basis of {toSplitEnum} : \n\n"


def lenAfterSplit_str(toSplitEnum, i, df_splitted):
    return f"   Length of DataSet after Splitting on basis on {toSplitEnum}{i + 1} is {len(df_splitted)}.\n"


def baseCond_str(boolvalue):
    if boolvalue == 2:
        return f"\n\nLength of DataSet is has become less than epsilon value of Length ({epsilonLen}). No need to Split Further!!"
    if boolvalue == 1:
        return f"\n\nEntropy for This DataSet is has become less than epsilon value of Entropy ({epsilonEntropy}). No need to Split Further!!"

    return f"\n\nWe have Reached to the bottom of the Decision Tree after Splitting all the Attributes.\n\n"


# function to calculate Significance of a split on basis of a variant of a particular Determining Attribute
def Significance_str(infoGain, n, N):
    sig = 0
    if N > 0:
        sig = infoGain * (n / N)
    return f"   Significance of This Split is : {sig}\n\n"


# function to calculate probability of Each Outcome in a particular dataset
def Probabilities_str(df):
    N = len(df)
    hashMap = dict(Counter(df["Output"]))

    info = ""
    for idx in possibleOutcomes:
        if idx in hashMap.keys():
            info += f"P({idx}) = {hashMap[idx] / N}\n"
        else:
            info += f"P({idx}) = 0.0\n"

    return info


if __name__ == "__main__":
    main()




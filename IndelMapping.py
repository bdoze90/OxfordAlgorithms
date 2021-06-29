"""Made for taking CLC Variant Output and representing it graphically"""


# Each item in the CLC dictionary (keyed on the gene/sequence ID) is a list of data with the following structure:
# Position
# Allele
# Reference Allele
# Coverage
# Frequency
# Ratio
CLCinfo = dict()


f = open("/Users/brianmendoza/Desktop/CLCinfoMatrix.txt")
for line in f:

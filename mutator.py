#!/usr/bin/python

#---------[ Info of this software - Begin ]---------#
# Author: Mehrad Mahmoudian
# Report bug: m.mahmoudian@gmail.com
#
# Version Number: 0.0.6.1
# Change Log:
#     - Add progressbar
#     - Add feature to check if the result file already
#        exists, and if the user is willing to overwrite
#        or specify new file name.
#     - Some comment modification
#     - TUI modification and uniformation (adding: Status, Caution, Warning and Message)
#
#---------[ Info of this software - End ]---------#


#--------[ Imports - Begin ]--------#
import argparse  # command line argument parser
import re # regular expression
from os import path # for checking existing directory
from os import makedirs # for creating directory
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
from math import ceil
import time
#--------[ Imports - End ]--------#

#---------[ Functions - Begin]----------#
def verifyEmail (inpt_str): # email validation function
  if len(inpt_str) > 7:
		if re.match("^.+\\@(\\[?)[a-zA-Z0-9\\-\\.]+\\.([a-zA-Z]{2,3}|[0-9]{1,3})(\\]?)$", inpt_str) != None:
			return True
		else:
			return False
	else:
		return False

def printError (inpt_field_name='Input data'): # prints errors related to the action
	print "\nWarning - The",inpt_field_name,"you entered was invalid, try again !\n"


def mkdir (dir_name):  # make directory if it does not exist in the given path
	if not path.exists(dir_name):
		makedirs(dir_name)

def writeFASTA (file_seq, file_name) :
	if len(file_seq) > 0 :
		output_handle = open(file_name, "w")
		SeqIO.write(file_seq, output_handle, "fasta")
		output_handle.close()
	else:
		print 'Internal Warning - No sequence is pushed in the writeFASTA function'

def getEntrezEmail () :
	correctEmailAddr = False
	while correctEmailAddr != True :
		EntrezEmail = raw_input("\nPlease enter a valid email address: { for example: sample@example.com }\n  > ")
		if(verifyEmail(EntrezEmail) == True):
			correctEmailAddr = True
		else:
			printError('Email address')
	return EntrezEmail

def getAccessionNumber () :
	correctAccessionNum = False
	while correctAccessionNum != True :
		AccessionNum = raw_input("\nPlease Type the Accession Number of the Protein: { for example: Q06187 }\n  > ")
		if len(AccessionNum)>2:
			correctAccessionNum = True
		else:
			printError('accession number')
	return AccessionNum

def getFastaFileName () :
	correctFileName = False
	while correctFileName != True :
		fastaFileName = raw_input("\nPlease Type the FASTA file name: { for example: results.fasta  or  result }\n  > ")
		if len(fastaFileName) > 0:  # if the input file name be not empty
			if fastaFileName[-6:len(fastaFileName)].lower() != '.fasta':  # if the provided file name does not contain .fasta extention, then add to it
				fastaFileName += '.fasta'
			correctFileName = True
		else:
			printError('FASTA file name')
	return fastaFileName

def checkFileExists (theFile) :
	# checks if the file exists, and if exists will return True
	try:
		with open(theFile): pass
		return True
	except IOError:
		return False

#---------[ Function - End ]----------#


#---------[ Main Code - Begin ]---------#
# [ getting argument from command line ]
parser = argparse.ArgumentParser(prog='Mutator')  # program name
parser = argparse.ArgumentParser(description='This is program is designed to create the mutation mentioned in genebank profile of a gene by getting the accession number and create a FASTA file output as a result.')  # program description
parser.add_argument('--version', action='version', version='%(prog)s 0.0.6.1')
parser.add_argument('--accession', action='store', help='Genebank accession number')
parser.add_argument('--email', action='store', help='Users email address for sending to Entrez')
parser.add_argument('--filename', action='store', help='Desired FASTA file name. FASTA extention is optional')

args = parser.parse_args()



# [ getting Accession Number from user ]
if args.accession is None :
	AccessionNum = getAccessionNumber ()
elif len(args.accession) < 2 :
	AccessionNum = getAccessionNumber ()
else:
	AccessionNum = args.accession


# [ getting email address from user ]
if args.email is None :
	Entrez.email = getEntrezEmail ()
elif verifyEmail(args.email) == False :
	Entrez.email = getEntrezEmail ()
else:
	Entrez.email = args.email

# [ getting FASTA file name from user ]
if args.filename is None :
	fastaFileName = getFastaFileName ()
elif len(args.filename) < 1 :
	fastaFileName = getFastaFileName ()
else:
	fastaFileName = args.filename
	if fastaFileName[-6:len(fastaFileName)].lower() != '.fasta':  # if the provided file name does not contain .fasta extention, then add to it
		fastaFileName += '.fasta'

fileOverWrite = ''
while fileOverWrite == '' :
	print "\nStatus  - Checking if the specified file exists:"
	if checkFileExists(fastaFileName) == False :  # check whether the file exists
		fileOverWrite = False
		print "            The file does not exist, we are good to go."
	else:
		fileOverWriteDecision = raw_input("Caution - The file you specified exists.\n            Do you want to specify a new file name? {type Y or N}\n  > ")
		if fileOverWriteDecision.lower() == "n" :
			fileOverWrite = True
		elif fileOverWriteDecision.lower() == "y":
			fastaFileName = getFastaFileName ()
		else:
			print "Your answer should be either Y or N."




# [ retrieving one record from Entrez service ]
print "\nStatus  - Retrieving data from Entrez service"
myErrorState = True
while myErrorState == True :
	handle = ""
	try:
		handle = Entrez.efetch(db="protein", rettype="gb", retmode="text", id=AccessionNum)
		myErrorState = False
	except Exception as inst:
		myErrorState = True
		print inst
		print "Warning - Seems the accession number is not correct, recheck the accession number and try again."
		AccessionNum = getAccessionNumber ()
if myErrorState == False :
	seq_record = SeqIO.read(handle, "genbank")  # use SeqIO.read when only one Seq
	handle.close()
else:
	sys.exit("""+-----------------------------------+
| Something unexpected happened !   |
|                                   |
| Please re-run the application     |
|    and if it happened again       |
|    please report it.              |
+-----------------------------------+""")


# [ extracting the sequence ]
print "\nStatus  - Extracting the variants and making mutations:"
masterSeq = seq_record.seq

myVariantList = []  # the variable for storing the mutant sequences and the related info for each
unknownRecords = []

featuresLength = len(seq_record.features)

for i in range(featuresLength):  # go through all features to find those who are variations
	myFeature = seq_record.features[i]
	if myFeature.type == 'Region':
		if 'region_name' in myFeature.qualifiers : # if the key (region_name) be found in the dictionary
			if myFeature.qualifiers['region_name'][0] == 'Variant' :  # if it was variation, then:
				myLocation = myFeature.location
				myNote = myFeature.qualifiers['note']
				if myNote[0][1:5] == ' -> ' :  # if it be a real variant and not something like 'Missing', then:
					if myNote[0][0] == masterSeq[myLocation.start] :  # if the note and sequence be exact match
						myVariantSeq = masterSeq[0:myLocation.start] + myNote[0][5] + masterSeq[myLocation.end:len(masterSeq)]  # inserting the mutated aminoacid in the original sequence
						myVariantTitle = "mutation: " + myNote[0][0:6] + " | location: " + str(myLocation.start) + ":" + str(myLocation.end)
						myVariant = SeqRecord(myVariantSeq, id=myNote[0][-10:-1], description=myVariantTitle)
						myVariantList.append(myVariant)  # appending the variant note and the mutated sequence in a list to wite them in a file later.
					else:
						unknownRecords.append('The location in notes does not match the aminoacid in the sequence ! \n the feature ID is ' + i)
				else:
					unknownRecords.append('The feature number' + str(i) + 'can not be processed.\nThe note of this record is: ' + myNote[0])
	# progressbar stuff
	myPercentageProgressBar = int(ceil(i*100/featuresLength))
	myProcessedProgressBar = int(ceil(myPercentageProgressBar*30/100))
	myRemainedProgressBar = 30 - myProcessedProgressBar
	time.sleep(0.01)
	sys.stdout.write("\r     [ " + myProcessedProgressBar*'#' + myRemainedProgressBar*' ' + ' ]   ' + str(myPercentageProgressBar) + "%")
	sys.stdout.flush()


# [ writing the mutated sequences in FASTA format in FASTA file ]
print "\nStatus  - Start writing Variants in Fasta file."
output_handle = open(fastaFileName, "w")
SeqIO.write(myVariantList, output_handle, "fasta")
output_handle.close()
print "\nStatus  - Writing Variants in Fasta file successfully completed."

# [ Showing unknown records in features ]
print "\nMessage - Below you can see the list of unprocessed qualifiers:"
for i in range(len(unknownRecords)) :
	print '          ->', unknownRecords[i]

#---------[ Main Code - End ]---------#

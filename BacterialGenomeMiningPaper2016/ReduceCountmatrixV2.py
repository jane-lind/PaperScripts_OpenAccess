#!/usr/bin/python2.7 -w

""" DESCRIPTION:
	The program reduces the count_matrices by a predefined average minimum count, 
	that has to be present within each row.
	0. Input argument and help program
	1. Printing overview of the arguments onto screen
	2. Print rows defining shared sets of genes.
 
INPUT FILE FORMATS:
	Input file: tab separated .txt file. Preferably BLAST output count matrix (See BLAST_profile.py).
  
Help:
	python ./reduce_countmatrix.py -h
  
Run example:
	python ReduceCountmatrixV2.py -i ./TestFiles/countMatrix_sorted_Test.txt -o ./ExtractionDir/
	
"""

# IMPORT MODULES
import os, datetime, getpass
from os.path import isfile, join
from sys import argv
from argparse import ArgumentParser
import re
import csv


# INPUT ARGUMENTS DESCRIPTION
# 0. Saving input arguments + help program
def ParseArguments():
	""" Handles program parameters and returns an argument class containing 
	all parameters """
	parser = ArgumentParser(description='Description of program: The program reduced the count_matrices by a predefined average minimum count. Each line is evaluated and only lines with an average count above the predefined minimum count will be generated.')
	parser.add_argument("-V", "--version", action="version", version="%(prog)s (version 0.1)")
	# Input parameters
	parser.add_argument("-i", type=str, dest="input_file", 
	default="", help="input matrix from BLASTprofile. Required input", required=True)
	parser.add_argument("-o", type=str, dest="output_dir", 
	default="", help="directory of your input file. Required input", required=True)
						
	args = parser.parse_args()
	
	# PRINT RUNNING DESCRIPTION 
	now = datetime.datetime.now()
	print '# ' + ' '.join(argv)
	print '# ' + now.strftime("%a %b %d %Y %H:%M")
	print '# USER: ' + getpass.getuser()
	print '# CWD: ' + os.getcwd()
	if os.name != "nt":
		print '# HOST: ' + ' '.join([ os.uname()[0] , os.uname()[4] , 
	       	                         os.uname()[2] , os.uname()[1] ])							
	return args



# ________SUB PROGRAMS___________
# 1. Printing overview of input arguments
def InputChecker(args):
	if os.path.os.path.isfile(args.input_file) == False:
		return sys.exit("\nThe input file do not exist: %s\nRerun the program with a new file entry.\n" %args.input_file)
	elif os.path.exists(args.output_dir) == False:
		return sys.exit("\nThe output directory do not exist: %s\nRerun the program with a new directory.\n" %args.output_dir)
	else:
		print '\n# COMMAND LINE PARAMETERS SET TO:'
		print "#\t Input file including path: ", args.input_file
		print "#\t Directory for saving output file: ", args.output_dir
		print ""



# 2A. Function for determining which row and column has the higher value in a given matrix (list of lists)
def FindMaximum(matrix, start_col):
	n = 0 # RowCounter
	maximum = int(0)
	for line in matrix:
		if max(line[start_col:]) > maximum:
			maximum = max(line[start_col:])
			column = line.index(maximum)
			row = n
		n += 1
	return(column, row, maximum)
	
# 2B. Compare a matrix (dimensions m x n) with a indexlist of length m. The function outputs a vector of length n, marking rows in the matrix with non-zero-entries where the vector has non-zero entries. The rows are given in reverse numerical order. the options start_col and start_row give the option of excluding rows and columns e.g. holding titles and labels
def FindRowsWithValues(matrix, indexlist, startcol=0, startrow=0):
	row_indexes = list()
	col = startcol
	for row in range(startrow, len(matrix)):
		for col in range(startcol, len(indexlist)):
			if indexlist[col] == int(0):
				pass
			elif indexlist[col] != int(0) and matrix[row][col] == int(0):
				pass
			elif indexlist[col] != int(0) and matrix[row][col] != int(0):
				row_indexes.append(int(row))
				break
			else:
				print "Something unforeseen has happened at row %(ROW)d and column %(COL)d!" % \
				{"ROW": row, "COL": col}
	row_indexes.sort()
	row_indexes.reverse()
	return(row_indexes)
	
# 2C: Count non-zero entries in a list
def CountNonzeroEntries(count_list, startcol=0):
	count = int(0)
	col = startcol
	for col in range(startcol, len(count_list)):
		if count_list[col] != 0:
			count += 1
	return(count)

# 2. Remove rows with an average count lower than given minimum
def RemoveExtraRows(args):
	f_matrix = open(args.input_file, "r")
	csv_matrix = list(csv.reader(f_matrix, delimiter='\t'))
	output_csv = open(args.output_dir+"/countMatrix_reduced.txt", 'w')
	f_matrix_reduced = csv.writer(output_csv, quoting=csv.QUOTE_ALL) #, dialect='excel-tab')
	
	f_matrix_reduced.writerow(csv_matrix[0]) # Write header to file
	del(csv_matrix[0]) # Remove the first line from the file
	# Convert csv_matrix to integers and strings
	for j in range(len(csv_matrix)):
		for i in range(1 , len(csv_matrix[j]) ):
			if csv_matrix[j][i] == '':
				csv_matrix[j][i] = int(0)
			else:
				csv_matrix[j][i] = int(csv_matrix[j][i])
	# Inititing loop
	n = 0
	species_counter = csv_matrix[0][1] # Counter for the species column
	print "Remaining rows to process:\n%s" % len(csv_matrix)	
	# Step 0: Make loop going through the lines one at a time
	while n < len(csv_matrix):
		# Convert all strings to intergers
		# for i in range(1 , len(csv_matrix[n]) ):
		# 	csv_matrix[n][i] = int(csv_matrix[n][i])
		# Check whether we are starting a new set.
		if csv_matrix[n][1] == species_counter:
			n += 1
		else:
			print "%s" % len(csv_matrix)
			# Reset variables
			group_matrix = csv_matrix[0:n]
			del(csv_matrix[0:n])
			n = 0
			# Run loop, delete max rows which don't have same number of entries as species-counter
			while len(group_matrix) > 0: # Write output to file if it passes QC and chew up the matrix as it goes along	
				max_data = FindMaximum(group_matrix, 2)
				species_hits_in_max_row = CountNonzeroEntries(group_matrix[max_data[1]], startcol=2)
				if species_hits_in_max_row == species_counter:
					f_matrix_reduced.writerow(group_matrix[max_data[1]]) # Write the maximum to the file
					deletionrow_indexes = FindRowsWithValues(group_matrix, group_matrix[max_data[1]], startcol=2, startrow=0)
					for del_index in deletionrow_indexes:
						del(group_matrix[del_index])
				else:
					del(group_matrix[max_data[1]])
				 
			# print CountNonzeroEntries(group_matrix[max_data[1]], startcol=2)
			species_counter = csv_matrix[0][1] # Reset species counter to new top row
			
	# Evaluating the last group (species_counter) of entries 
	group_matrix = csv_matrix
	while len(group_matrix) > 0: # Write output to file if it passes QC and chew up the matrix as it goes along	
		max_data = FindMaximum(group_matrix, 2)
		species_hits_in_max_row = CountNonzeroEntries(group_matrix[max_data[1]], startcol=2)
		if species_hits_in_max_row == species_counter:
			f_matrix_reduced.writerow(group_matrix[max_data[1]]) # Write the maximum to the file
			deletionrow_indexes = FindRowsWithValues(group_matrix, group_matrix[max_data[1]], startcol=2, startrow=0)
			for del_index in deletionrow_indexes:
				del(group_matrix[del_index])
		else:
			del(group_matrix[max_data[1]])
	
	# Remember to end with evaluating the last term
	f_matrix.close()
	output_csv.close()
        
    


							# 


# ______ MAIN PROGRAM______
def main(argv):
	args = ParseArguments()																	# 0.
	InputChecker(args)																				# 1.
	RemoveExtraRows(args)																				# 2.



if __name__ == '__main__':
	main(argv)
	print "\n### Your program has finished! ###\n"
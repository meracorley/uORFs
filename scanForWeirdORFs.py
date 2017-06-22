def main():
	#Script to tell me where ORFs starting with a non-AUG start codon are in a given sequence.
	#It takes a sequence, by itself or in fasta format, and gives the relative position of ORFs in the sequence.
	#Doesn't distinguish between uORFs and main CDS because you can't without other gene annotation info.
	#Big assumption: the kozak sequence strength surrounding a non-AUG start is its AUG start Kozak strength * eff. of non-AUG start
	#Limitations: Input sequence has to all be on one line.
	#TO RUN: python scanForORFs.py seq.txt/.fa

	import sys, os
	storeCodons = [] #Reticulocyte translation initiation efficiencies relative to ATG. Taken from Peabody 1989 
	weirdStartCodons = {"ATG":1., "ACG":.84, "GTG":.36, "TTG":.39, "CTG":.82, "AGG":.17, "AAG":.14, "ATA":.59, "ATC":.47, "ATT":.67}
	kozakStrengthDict = {}
	nodererFile = open("/home/mcorley/kozakStrengths.txt", 'r') #Kozak strengths from Noderer et al
	noderer = nodererFile.read().splitlines()
	for line in noderer: #Putting strength of each kozak seq into a dictionary
		data = line.split()
		kozak = data[0]
		strength = data[1]
		kozakStrengthDict[kozak.replace('U','T')] = strength
	nodererFile.close()	

	if os.path.isfile(sys.argv[1]): #can input a file or a sequence	
        	infile = open(sys.argv[1], 'r')
        	seq = infile.read().splitlines()
		if seq[0][0]==">":
			seq = seq[1].strip('\n').upper()
			seq = seq.replace('U','T') #I use Ts instead of U. Personal preference.
		else:
			seq=seq[0].strip('\n').upper()
			seq = seq.replace('U','T')
		infile.close()
	else:
		seq = sys.argv[1].strip('\n').upper()
		seq = seq.replace('U','T')
	i = 0
	cdsStart = len(seq)-2
	if len(sys.argv)>2 and is.integer(int(sys.argv[2])): #can give CDS start at the second argument
		cdsStart = int(sys.argv[2])
	while i < cdsStart:
		codon = seq[i:i+3] 
		if codon in weirdStartCodons:
			#print "AUG at", i+1
			noStopCodon = True
			j = i+3
			while j < len(seq)-2:
				nextcodon = seq[j:j+3]
				storeCodons.append(nextcodon)
				if nextcodon=="TGA" or nextcodon=="TAA" or nextcodon=="TAG":
					noStopCodon = False
					#print storeCodons #printing the uORF
					#print seq[j-len(storeCodons)*3-5:j-len(storeCodons)*3+6] #printing the Kozak seq
					#Estimate the Kozak strength of the non-ATG uORF
					kozak = seq[i-6:i] + "ATG" + seq[i+3:i+5] #the kozakseq if start codon were an ATG
					kozakStrength = "Incomplete"
					if len(kozak)==11: #some initiation sites start too close to end to have a full kozak
						kozakStrength = weirdStartCodons[codon] * float(kozakStrengthDict[kozak])
					print i+1 ,j+3, kozakStrength #returns the beginning pos, end pos and Kozak seq strength of each ORF
					storeCodons = []
					#print seq[i:i+3], seq[j:j+3], ((j+3)-i)%3 
					break
				j+=3

			#if noStopCodon:
			#	print "AUG at", i+1
			i+=3		
		else:
			i+=1

main()

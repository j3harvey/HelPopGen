import argparse
parser = argparse.ArgumentParser()
parser.add_argument("threshold",
                    help="Minor allele frequency cutoff as a string",
                    nargs='?',
                    type=string,
                    default=0.0)
args = parser.parse_args()
THRES = args.threshold

def unique(I):
	"""
	Yields unique elements of iterable I in the order they first appear
	"""
	i = iter(I)
	y = []
	while 1:
		x = i.next()
		if x not in y:
			yield x
			y.append(x)


with open("../../data/HmelDmel_OrthologyGroups_FINAL_231013.csv", 'r') as f:
	headers = f.readline().rstrip().split('\t')
	immunity_genes = [dict(zip(headers, l.rstrip().split('\t'))) for l in f.readlines()]
	orthology_groups = list(unique([g['Orthology Group'] for g in genes]))
	groups_map = {gp : n+1 for (n, gp) in enumerate(orthology_groups)}
	genes_map = {g['GeneID'] : g['Orthology Group'] for g in genes}
	hmel_immunity_genes = [g for g in immunity_genes if g['GeneID'].startswith('HMEL')]
	dmel_immunity_genes = [g for g in immunity_genes if g['GeneID'].startswith('FBgn')]


with open("../../results/paper" + THRES + ".csv", 'r') as f:
	headers = f.readline().rstrip().split('\t')
	hmel_genes = [dict(zip(headers, l.rstrip().split('\t'))) for l in f.readlines()]


with open("../../data/s027_all.csv", 'r') as f:
	headers = f.readline().rstrip().split('\t')
	dmel_genes = [dict(zip(headers, l.rstrip().split('\t'))) for l in f.readlines()]

# csv writer for rows of the MKtest input file
mkin_headers = ["D_N",
		"Ln",
		"in_P_N",
		"Ln",
		"D_S",
		"Ls",
		"P_S",
		"Ls",
		"ingroup_sequences",
		"chr",
		"class",
		]
mkin_writer = DictWriter(open("mkin_t_immunity.csv", 'w'), headers)

for gene in hmel_genes:
	# write MK table giving genes a class of groups_map(group)
	gene["chr"] = 0
	# classes of 1 or more for homology groups, otherwise class is 0
	gene['class'] = genes_map.get(gene[''], 0)
	# calculate numbers of sites
	gene['Ln'] = float(gene['NSsites'])/3
	gene['Ls'] = float(gene['Ssites'])/3
	row = {k: gene(k) for k in mkin_headers}
	if row
	mkin_writer.writerow(row)

###
#
# Run MKtest
#

subprocess.call(['../../MKtest-2.0/MKtest', '-f2', 'mktest_immunity' + THRES + '.csv', '-o', 'immunity' + THRES + '.mkout'])



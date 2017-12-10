import sys
import argparse
import pkg_resources
gffutils_version = pkg_resources.get_distribution("gffutils").version
if float(gffutils_version) >= 0.899:
	import gffutils
else:
	print("This script requires gffutils version >= 0.9")
	exit()

parser = argparse.ArgumentParser()
parser.add_argument("file",help="The path to the GFF3 file to convert")
args = parser.parse_args()

db_path=args.file+".gffutils.db"
sys.stderr.write("Reading GFF3 file: "+args.file+"\n")
sys.stderr.write("Coverting to gffutils sqlite database: "+db_path+"\n")
sys.stderr.flush()
db = gffutils.create_db(args.file, db_path, force=True,merge_strategy="create_unique")
sys.stderr.write("Done converting. Now printing modified GFF3 to stdout...\n")
sys.stderr.flush()

print("##gff-version 3")
for g in db.features_of_type("gene"):
	##print g ##Print the gene line
	##Nevermind. Can't just print the gene line, as there may be multiple mRNAs, but one gene, which EVM can't handle
	
	mRNAs = db.children(g.id,featuretype='mRNA')
	for m in mRNAs:
		##Have to fake a gene line per mRNA
		gene_id = m.id+"-gene"
		gene_attributes = "ID="+gene_id+";Name="+str(m.id)
		gene_string = '\t'.join([m.chrom,m.source,"gene",str(m.start),str(m.stop),m.score,m.strand,".",gene_attributes])
		print(gene_string)
		## Done faking gene line
	
		m.attributes["ID"] = m.id ##<< Add a unique ID to the mRNA. Needed for EVM.  GFFUtils automatically makes it unique.
		m.attributes["Parent"] = gene_id ##Link it to the parent fake gene
		print m ##Print the mRNA line
		exons = db.children(m.id,featuretype='exon')	
		for e in exons:
			if len(e.attributes['Parent']) == 1:
				e.attributes["ID"] = e.id ##<< Add a unique ID to the exon. Needed for EVM.
				print e ##Print the exon line
			if len(e.attributes['Parent']) > 1:
				##Have to fake the id, to make it an independent exon
				ID = e.id+'_of_'+m.id ##Print the faked exon line
				new_attributes = "ID="+ID
				for a_key in e.attributes.keys():
					if a_key == "Parent":
						new_attributes = new_attributes+";Parent="+m.id
					else:
						new_attributes = new_attributes+";"+a_key+"="+str(e.attributes[a_key][0])
				feature_string = '\t'.join([e.chrom,e.source,e.featuretype,str(e.start),str(e.stop),e.score,e.strand,".",new_attributes])
				print(feature_string)

		CDSs = db.children(m.id,featuretype='CDS')	
		for c in CDSs:
			if len(c.attributes['Parent']) == 1:
				c.attributes["ID"] = c.id ##<< Add a unique ID to the exon. Needed for EVM.
				print c ##Print the CDS line
			if len(c.attributes['Parent']) > 1:
				##Have to fake the id, to make it an independent CDS
				ID = c.id+'_of_'+m.id ##Print the faked CDS line
				new_attributes = "ID="+ID
				for a_key in c.attributes.keys():
					if a_key == "Parent":
						new_attributes = new_attributes+";Parent="+m.id
					else:
						new_attributes = new_attributes+";"+a_key+"="+str(c.attributes[a_key][0])
				feature_string = '\t'.join([c.chrom,c.source,c.featuretype,str(c.start),str(c.stop),c.score,c.strand,c.phase,new_attributes])
				print(feature_string)

sys.stderr.write("Conversion complete.\n")

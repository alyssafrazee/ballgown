## development tests for ballgown package
## alyssa frazee 25 sep 2013

MAC = R.Version()$os == 'darwin10.8.0'
CLUSTER = R.Version()$os == 'linux-gnu'

workdir = ifelse(MAC, 
	"~/Google Drive/hopkins/research/_ballgown/testlogs",
	"/home/bst/student/afrazee/ballgown/dev")

# set up file for log of tests
test_date = paste(strsplit(date(), split=' ')[[1]][2:5], collapse='_')
sink(paste0(workdir, '/testlog_',test_date))

# install the package
library(devtools)
cat(date())
cat(': installing package...\n')
cat(try(install_github("ballgown", "alyssafrazee", subdir="ballgown")))
cat('\n\n')
library(ballgown)

# test constructor function
dataDir = ifelse(MAC, 
	'~/Documents/hopkins/research/_RNAseq-simulation',
	'/amber2/scratch/jleek/orbFrontal/scripts/simulations/ballgown')
samplePattern = 'sample'
files = list.files(dataDir)
pData = data.frame(
	dirname = files[which(substr(files,1,2) == "sa")],
	group = c("A",rep("B",9),"A","B",rep("A",8)))
cat(date())
cat(': testing constructor function...\n')
print(try(simgown <- ballgown(dataDir = dataDir, samplePattern = samplePattern, pData = pData)))
cat('\n\n')

# test subset method
cat(date())
cat(': testing subset method...\n')
print(try(simgown <- subset(simgown, chr==22)))
cat('\n\n')


sink()


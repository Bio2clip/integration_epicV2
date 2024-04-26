library(readr)

# Importing SampleSheet from CSV file
epicv2 <- read_delim(
  file = "SAMPLES.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Importing CNV data
load.data <- conumee2.0::CNV.import(
  array_type = "EPICv2",
  directory = "DIR",
  sample_sheet = epicv2
)

# Loading exclusion and detail regions
data(exclude_regions)
data(detail_regions)

# Selecting detail regions
OLD_Region <- detail_regions
myGR <- GRanges(
  seqnames = c("chr1"),
  ranges = IRanges(start = 241660903, end = 241683055),
  strand = "-",
  name = 'FH',
  thick = IRanges(start = 241660903, end = 241683055)
)
NEW_REGION <- c(OLD_Region, myGR)

# Creating annotation
anno <- CNV.create_anno(
  array_type = "EPICv2",
  exclude_regions = exclude_regions,
  detail_regions = NEW_REGION
)

# Loading CNV data
load.data <- CNV.load(load.data)

x <- CNV.fit(query = load.data[SAMPLE],ref = load.data[SAMPLE_CONTROL], anno)
x <- CNV.segment(x)
x <- CNV.detail(x)
CNV.genomeplot(x)

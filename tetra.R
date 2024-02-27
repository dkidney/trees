# install.packages('tidyverse')
# install.packages('BiocManager')
# BiocManager::install("ggtree")
# install.packages('phytools')
# install.packages('deeptime')
# install.packages('TreeTools')

suppressPackageStartupMessages({
	library(magrittr)
	library(tidyverse)
	library(ggtree)
	library(phytools)
})

# browseVignettes("ggtree")
# https://yulab-smu.top/treedata-book/
# https://4va.github.io/biodatasci/r-ggtree.html
# https://cran.r-project.org/web/packages/deeptime/vignettes/phylogenies.html
# https://github.com/willgearty/deeptime
# https://bioinformatics.stackexchange.com/questions/20677/plot-phylogenetic-tree-from-list-of-edges

tetra = ape::read.tree('tetra.tree')
tetra
tetra %>% as_tibble %>% print(n=Inf)
tetra %>% as_tibble %>% filter(!node %in% parent)
tetra %>% as_tibble %>% filter(node %in% parent)

methods(class=class(tetra))
Nnode(tetra)
Ntip(tetra)
Nedge(tetra)
degree(tetra)
ltt(tetra)
tetra %>% ggtree()
plt = tetra %>%
	# ggtree(aes(color=group)) +
	ggtree() +
	geom_vline(aes(xintercept=-358.9), lty=2, col='grey75') +
	geom_vline(aes(xintercept=-298.9), lty=2, col='grey75') +
	geom_vline(aes(xintercept=-252), lty=2, col='grey75') +
	geom_vline(aes(xintercept=-201.4), lty=2, col='grey75') +
	geom_vline(aes(xintercept=-145), lty=2, col='grey75') +
	geom_vline(aes(xintercept=-66), lty=2, col='grey75') +
	geom_vline(aes(xintercept=-23.03), lty=2, col='grey75') +
	geom_vline(aes(xintercept=-2.58), lty=2, col='grey75') +
	geom_vline(aes(xintercept=-0), lty=2, col='grey75') +
	# theme_tree2() +
	# geom_text(aes(label=label), hjust=-.1) +
	geom_tiplab(hjust = 0, geom='label', size=3) +
	geom_nodelab(hjust = 1.3, geom='label', size=3) +
	geom_nodepoint() +
	geom_tippoint() +
	# geom_hilight(node=33, fill="gold") +
	# geom_label(aes(x=branch, label=label)) +
	# geom_label(aes(label=label2)) +
	deeptime::coord_geo(
		xlim = c(-358.9, 50),
		ylim = c(0, ape::Ntip(tetra) + 1),
		dat = list('stages', "epochs", "periods"),
		pos = list("bottom", "bottom", "bottom"),
		size = list(2, 3, 5),
		neg = TRUE,
		abbrv = TRUE,
		center_end_labels = TRUE
	) +
	scale_x_continuous(breaks = seq(-1000, 0, 50),
					   labels = -seq(-1000, 0, 50),
					   expand = c(1000, 100)) +
	theme_tree2()
# print(plt)
plt = revts(plt)
print(plt)









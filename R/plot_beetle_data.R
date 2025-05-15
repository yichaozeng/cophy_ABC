# data needed to visualize the beetle cophylogeny
ass_mat_trans <- NULL
ass_mat <- cophy_real[[1]][[3]]
for(i in 1:nrow(ass_mat)){
  for (j in 1:ncol(ass_mat)){
    if (ass_mat[i,j] == 1)
      ass_mat_trans <- rbind(ass_mat_trans, c(rownames(ass_mat)[i], colnames(ass_mat)[j]))
  }
}

tree1 <- cophy_real[[1]][[1]]
tree2 <- cophy_real[[1]][[2]]
net <- cophy_real[[1]][[3]]

# we use ape / phytools to visualize the associations, as it is the only package that properly displays the links
cophylo_plot <- cophylo(tree1, tree2, assoc = ass_mat_trans, rotate = F)

pdf(paste('cophy_ABC/R/real_data/links.pdf', sep = ''), width = 20, height = 20, pointsize = 12)
plot.cophylo(cophylo_plot,link.lwd=1,
             link.lty="solid",link.col=make.transparent("black",1))
dev.off()

# we use the treeducken package to visualize the two trees, as it is the only package that can properly display the root edge of a phylogeny
cophylo_plot_treeducken <- list(host_tree=tree1, symb_tree=tree2, association_mat=net)
class(cophylo_plot_treeducken) <- "cophy"

pdf(paste('cophy_ABC/R/real_data/trees_with_roots.pdf', sep = ''), width = 40, height = 20, pointsize = 12)
plot.cophy(cophylo_plot_treeducken)
dev.off()

# we use the ape package to visualize the scale bar
pdf(paste('cophy_ABC/R/real_data/scale_bar.pdf', sep = ''), width = 10, height = 20, pointsize = 12)
plot.phylo(tree1, root.edge = T, show.tip.label = F)
axis(side=1)
dev.off()

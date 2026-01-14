# 加载包
library(ape)
library(ggtree)

# 读取 Newick 树
newick_str <- "(Apis_mellifera:0.361049,(((Gryllus_bimaculatus:0.343318,Blattella_germanica:0.34685):0.0469324,Acyrthosiphon_pisum:0.582244):0.0217689,((Bombyx_mori:0.242494,Helicoverpa_armigera:0.359651):0.128642,(Tribolium_castaneum:0.322842,(Drosophila_mojavensis:0.451666,Aedes_aegypti:0.35007):0.110205):0.0170747):0.0114513):0.0352473);"
tree <- read.tree(text = newick_str)

# 绘制从上往下的竖排树
p <- ggtree(tree,
            layout        = "rectangular",   # 矩形树
            # branch.length = "none",          # 若要真实分支长度，可移除此行
            direction     = "down") +        # 让树从上向下延伸
  geom_tiplab(
    hjust = 0, 
    size  = 3
  ) +
  theme_tree2() +
  scale_x_continuous(name = "Million years ago") +
  xlab(NULL) + ylab(NULL)

print(p)

ggsave("tree.pdf",p ,height = 6, width=1,dpi=300)
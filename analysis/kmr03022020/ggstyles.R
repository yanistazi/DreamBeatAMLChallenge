### ggplot styles
library(ggplot2)
theme1 = theme_bw() + theme(text=element_text(face="bold", size=20))
theme0 = theme_bw() + theme(text=element_text(face="bold", size=10))
noleg = theme(legend.position = "none")
topleg = theme(legend.position = "top")
angle45 = theme(axis.text.x = element_text(angle = 45, hjust = 1))
noxtitle = theme(axis.title.x = element_blank())
noytitle = theme(axis.title.y = element_blank())
noxlabel = theme(axis.text.x = element_blank())
noylabel = theme(axis.text.y = element_blank())
nolegtitle = theme(legend.title = element_blank())

legtwoshape = guides(shape=guide_legend(ncol=2))
legtwofill = guides(fill=guide_legend(ncol=2))


bold15  = theme(axis.text = element_text(size=15,face="bold"))
bold20  = theme(axis.text = element_text(size=20,face="bold"))


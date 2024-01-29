library(tidyverse)

N_fmu = 25
dV = (10/N_fmu)^2
x = seq(0,10, length.out=51)
y = seq(0,10, length.out=51)

u_iebt = read.delim("~/codes/libpspm/iebt2d_u.txt")
u1_iebt = read.delim("~/codes/libpspm/iebt2d_u1.txt")
dx_iebt = read.delim("~/codes/libpspm/iebt_wave_x.txt", header=F)
dy_iebt = read.delim("~/codes/libpspm/iebt_wave_y.txt", header=F)

u_ifmu = read.delim("~/codes/libpspm/ifmu2d_u.txt")
u1_ifmu = read.delim("~/codes/libpspm/ifmu2d_u1.txt")
dx_ifmu = read.delim("~/codes/libpspm/ifmu_wave_x.txt", header=F)
dy_ifmu = read.delim("~/codes/libpspm/ifmu_wave_y.txt", header=F)
plot_model(u_ifmu, u1_ifmu, dx_ifmu, dy_ifmu, dV, title="IFMU")

u_abm = read.delim("~/codes/libpspm/abm2d_u.txt")
u1_abm = read.delim("~/codes/libpspm/abm2d_u1.txt")
dx_abm = read.delim("~/codes/libpspm/abm_wave_x.txt", header=F)
dy_abm = read.delim("~/codes/libpspm/abm_wave_y.txt", header=F)

plot_model = function(u, u1, dx, dy, dV=1, title){
  pp0 = u %>% ggplot(aes(x=x.0., y=x.1., fill=u*dV)) + 
    geom_rect(aes(xmin=0, xmax=10, ymin=0, ymax=10), fill="#132B43")+
    geom_tile(width=0.4, height=0.4) + 
    geom_point(aes(x=2,y=4), col="red") + 
    xlim(c(0,10)) + ylim(c(0,10))
  
  pp1 = u1 %>% ggplot(aes(x=x.0., y=x.1., fill=u*dV)) + 
    geom_rect(aes(xmin=0, xmax=10, ymin=0, ymax=10), fill="#132B43")+
    geom_tile(width=0.4, height=0.4) + 
    geom_point(aes(x=2,y=4), col="red") + 
    geom_point(aes(x=2+1*1,y=4+2*1), col="yellow") + 
    xlim(c(0,10)) + ylim(c(0,10)) 
  
  p1 = autoplot(zoo::zoo(t(dx), order.by = seq(0,10, length.out=51)), facet = NULL) + 
    scale_color_manual(values = 
                         scales::seq_gradient_pal("pink","#08306b")(seq(0,1,length.out=nrow(dx)))) + 
    # scale_y_reverse()+
    guides(colour="none") 
  p2 = autoplot(zoo::zoo(t(dy), order.by = seq(0,10, length.out=51)), facet = NULL) + 
    scale_color_manual(values = 
                         scales::seq_gradient_pal("pink","#08306b")(seq(0,1,length.out=nrow(dx)))) + 
    # scale_color_brewer("seq") + 
    # scale_y+reverse()+
    guides(colour="none")
  
  # print(
  # cowplot::plot_grid(p1+theme_classic()+labs(x=""),
  #                    pp0+guides(fill="none")+theme_classic()+labs(x="t=0", y=""),
  #                    pp1+guides(fill="none")+theme_classic()+labs(x="t=1", y=""),
  #                    p2+theme_classic()+labs(x=""), 
  #                    align="hv", axis = "lbrt")
  # )

  # cairo_pdf(paste0(title, "_2d_wave.pdf"))
  print(
    cowplot::plot_grid(
      p2+theme_classic()+labs(x="", y="U(y)")+scale_y_reverse()+coord_flip()+theme(axis.text.y = element_blank())+scale_x_continuous(position="top")+ggtitle(paste0("\n",title)), 
      pp1+guides(fill="none")+theme_classic()+labs(x="x", y="y")+ggtitle("u(x,y, t=1)"),
      pp0+guides(fill="none")+theme_classic()+labs(x="", y="")+ggtitle("u(x,y, t=0)"),
      p1+theme_classic()+labs(x="", y="U(x)")+scale_y_reverse()+theme(axis.text.x = element_blank())+scale_x_continuous(position="top"),
      align="hv", axis = "lbrt", greedy = T)
  )
  # dev.off()
  list(px=p1+theme_classic()+labs(x="x",y="U(x)")+geom_vline(xintercept=2+1*1, col="yellow2")+ylim(c(0,0.5)), 
       py=p2+theme_classic()+labs(x="y",y="U(y)")+geom_vline(xintercept=4+2*1, col="yellow2")+ylim(c(0,0.5)),
       pp0=pp0+guides(fill="none")+theme_classic()+labs(x="x", y="y"), 
       pp1=pp1+guides(fill="none")+theme_classic()+labs(x="x", y="y"))
}

l_iebt = plot_model(u_iebt, u1_iebt, dx_iebt, dy_iebt, title="IEBT")
l_ifmu = plot_model(u_ifmu, u1_ifmu, dx_ifmu, dy_ifmu, dV, title="IFMU")
l_abm = plot_model(u_abm, u1_abm, dx_abm, dy_abm, title="ABM")

cairo_pdf("combined_2d_wave.pdf")
cowplot::plot_grid(
  l_iebt$px+ggtitle("IEBT"), l_iebt$py, l_iebt$pp1,
  l_ifmu$px+ggtitle("IFMU"), l_ifmu$py, l_ifmu$pp1,
  l_abm$px+ggtitle("ABM"), l_abm$py, l_abm$pp1,
  align="hv", byrow = F
)
dev.off()

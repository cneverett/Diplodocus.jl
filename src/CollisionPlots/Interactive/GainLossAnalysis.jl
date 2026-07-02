function InteractiveBinaryGainLossPlot(Output::Tuple)

    GLMakie.activate!(inline=false)

    (name1,name2,name3,name4,mu1,mu2,mu3,mu4,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,p4_low,p4_up,p4_grid,p4_num,u4_grid,u4_num,h4_grid,h4_num) = Output[1]
    
    p1_r = DC.bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = DC.deltaVector(p1_r);
    p1_m = DC.meanVector(p1_r);
    p2_r = DC.bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = DC.deltaVector(p2_r);
    p2_m = DC.meanVector(p2_r);
    p3_r = DC.bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_d = DC.deltaVector(p3_r);
    p3_m = DC.meanVector(p3_r);
    p4_r = DC.bounds(p4_low,p4_up,p4_num,p4_grid);
    p4_d = DC.deltaVector(p4_r);
    p4_m = DC.meanVector(p4_r);
    
    u1_r = DC.bounds(-1.0,1.0,u1_num,u1_grid);
    u2_r = DC.bounds(-1.0,1.0,u2_num,u2_grid);
    u3_r = DC.bounds(-1.0,1.0,u3_num,u3_grid);
    u4_r = DC.bounds(-1.0,1.0,u4_num,u4_grid);
    u1_d = DC.deltaVector(u1_r);
    u2_d = DC.deltaVector(u2_r);
    u3_d = DC.deltaVector(u3_r);
    u4_d = DC.deltaVector(u4_r);
    u1_m = DC.meanVector(u1_r);
    u2_m = DC.meanVector(u2_r);
    u3_m = DC.meanVector(u3_r);
    u4_m = DC.meanVector(u4_r);
    h1_r = DC.bounds(0.0,2.0,h1_num,h1_grid);
    h2_r = DC.bounds(0.0,2.0,h2_num,h2_grid);
    h3_r = DC.bounds(0.0,2.0,h3_num,h3_grid);
    h4_r = DC.bounds(0.0,2.0,h4_num,h4_grid);
    h1_d = DC.deltaVector(h1_r);
    h2_d = DC.deltaVector(h2_r);
    h3_d = DC.deltaVector(h3_r);
    h4_d = DC.deltaVector(h4_r);
    h1_m = DC.meanVector(h1_r);
    h2_m = DC.meanVector(h2_r);
    h3_m = DC.meanVector(h3_r);
    h4_m = DC.meanVector(h4_r);

    fig = Figure(size = (1000,600))
    
    gl1 = GridLayout(fig[3,1:4],tellwidth = false)
    gl2 = GridLayout(fig[1:2,3:4],tellwidth = false)
    subgl1 = GridLayout(gl1[1,1],tellwidth = false)
    subgl2 = GridLayout(gl2[1,1],tellwidth = false)

    sg1 = SliderGrid(subgl1[1,1],
    (label = "p1", range = 1:length(p1_m), startvalue = 1, update_while_dragging = false),
    (label = "u1", range = 1:length(u1_m), startvalue = 1, update_while_dragging = false),
    (label = "h1", range = 1:length(h1_m), startvalue = 1, update_while_dragging = false)
    )
    sg2 = SliderGrid(subgl1[1,2],
    (label = "p2", range = 1:length(p2_m), startvalue = 1, update_while_dragging = false),
    (label = "u2", range = 1:length(u2_m), startvalue = 1, update_while_dragging = false),
    (label = "h2", range = 1:length(h2_m), startvalue = 1, update_while_dragging = false)
    )
    sg3 = SliderGrid(subgl1[1,3],
    (label = "u3", range = 1:length(u3_m), startvalue = 1, update_while_dragging = false),
    (label = "h3", range = 1:length(h3_m), startvalue = 1, update_while_dragging = false)
    )
    sg4 = SliderGrid(subgl1[1,4],
    (label = "u4", range = 1:length(u4_m), startvalue = 1, update_while_dragging = false),
    (label = "h4", range = 1:length(h4_m), startvalue = 1, update_while_dragging = false)
    )

    cbAngAvg = Checkbox(subgl2[2,1],checked=false)
    cbAnalytic = Checkbox(subgl2[3,1],checked=false)
    Label(subgl2[2,2],"Angle Averaged",halign=:left)
    Label(subgl2[3,2],"Analytic ISO Spectrum",halign=:left)

    p1_idx = sg1.sliders[1].value
    u1_idx = sg1.sliders[2].value
    h1_idx = sg1.sliders[3].value
    p2_idx = sg2.sliders[1].value
    u2_idx = sg2.sliders[2].value
    h2_idx = sg2.sliders[3].value
    u3_idx = sg3.sliders[1].value
    h3_idx = sg3.sliders[2].value
    u4_idx = sg4.sliders[1].value
    h4_idx = sg4.sliders[2].value

    p1_val = @lift(log10.(p1_m[$p1_idx]))
    u1_val = @lift(u1_m[$u1_idx])
    h1_val = @lift(h1_m[$h1_idx])
    p2_val = @lift(log10.(p2_m[$p2_idx]))
    u2_val = @lift(u2_m[$u2_idx])
    h2_val = @lift(h2_m[$h2_idx])
    Analytic = @lift($(cbAnalytic.checked) == true)

    label = @lift("p1: [$(round(p1_r[$p1_idx], sigdigits=3)),$(round(p1_r[$p1_idx+1], sigdigits=3))], u1: [$(round(u1_r[$u1_idx],sigdigits=2)),$(round(u1_r[$u1_idx+1],sigdigits=2))], h1: [$(round(h1_r[$h1_idx],sigdigits=2)),$(round(h1_r[$h1_idx+1],sigdigits=2))] 
    \n 
    p2: [$(round(p2_r[$p2_idx], sigdigits=3)),$(round(p2_r[$p2_idx+1], sigdigits=3))], u2: [$(round(u2_r[$u2_idx],sigdigits=2)),$(round(u2_r[$u2_idx+1],sigdigits=2))], h2: [$(round(h2_r[$h2_idx],sigdigits=2)),$(round(h2_r[$h2_idx+1],sigdigits=2))]
    \n
    u3: [$(round(u3_r[$u3_idx],sigdigits=2)),$(round(u3_r[$u3_idx+1],sigdigits=2))], h3: [$(round(h3_r[$h3_idx],sigdigits=2)),$(round(h3_r[$h3_idx+1],sigdigits=2))]
    \n
    u4: [$(round(u4_r[$u4_idx],sigdigits=2)),$(round(u4_r[$u4_idx+1],sigdigits=2))], h4: [$(round(h4_r[$h4_idx],sigdigits=2)),$(round(h4_r[$h4_idx+1],sigdigits=2))] 
    ") 

    Label(subgl2[4,1:2], label)

    GainMatrix3 = Output[2]
    GainMatrix3Avg = dropdims(sum(GainMatrix3,dims=(2,3,5,6,8,9)),dims=(2,3,5,6,8,9)) #= / sum(u3_d) * u3_d[1] =# / sum(u1_d) * u1_d[1] / sum(u2_d) * u2_d[1]# / (u1_num*h1_num*u2_num*h2_num*u3_num*h3_num)
    GainMatrix4 = name1==name2 ? Output[2] : Output[3]
    GainMatrix4Avg = dropdims(sum(GainMatrix4,dims=(2,3,5,6,8,9)),dims=(2,3,5,6,8,9)) #= / sum(u4_d) * u4_d[1] =# / sum(u1_d) * u1_d[1] / sum(u2_d) * u2_d[1]#/ (u1_num*h1_num*u2_num*h2_num*u4_num*h4_num)
    LossMatrix1 = Output[4]
    LossMatrix1Avg = dropdims(sum(LossMatrix1,dims=(2,3,5,6)),dims=(2,3,5,6)) #= / sum(u1_d) * u1_d[1] =# / sum(u2_d) * u2_d[1]#/ (u1_num*h1_num*u2_num*h2_num)
    LossMatrix2 = name1==name2 ? Output[4] : Output[5]
    LossMatrix2Avg = dropdims(sum(LossMatrix2,dims=(2,3,5,6)),dims=(2,3,5,6)) / sum(u1_d) * u1_d[1] #= / sum(u2_d) * u2_d[1] =##/ (u1_num*h1_num*u2_num*h2_num)
    #GainMatrix3 = Output[2]
    #GainMatrix4 = Output[3]
    #LossMatrix1 = Output[4]
    #LossMatrix2 = Output[5]
    #GainMatrix3 = dropdims(sum(Output[2],dims=(2,3,5,6,8,9)),dims=(2,3,5,6,8,9));
    #GainMatrix4 = dropdims(sum(Output[3],dims=(2,3,5,6,8,9)),dims=(2,3,5,6,8,9));
    #LossMatrix1 = dropdims(sum(Output[4],dims=(2,3,5,6)),dims=(2,3,5,6));
    #LossMatrix2 = dropdims(sum(Output[5],dims=(2,3,5,6)),dims=(2,3,5,6));


    # Loss Matrix   
    p1_loss = @lift($(cbAngAvg.checked) ? log10.(LossMatrix1Avg[$p1_idx,$p2_idx]) : log10.(LossMatrix1[$p1_idx,$u1_idx,$h1_idx,$p2_idx,$u2_idx,$h2_idx]))
    p2_loss = @lift($(cbAngAvg.checked) ? log10.(LossMatrix2Avg[$p2_idx,$p1_idx]) : log10.(LossMatrix2[$p2_idx,$u2_idx,$h2_idx,$p1_idx,$u1_idx,$h1_idx]))
    # Gain Matrix
    p3_gain = @lift($(cbAngAvg.checked) ? log10.(GainMatrix3Avg[2:end-1,$p1_idx,$p2_idx]) : log10.(GainMatrix3[2:end-1,$u3_idx,$h3_idx,$p1_idx,$u1_idx,$h1_idx,$p2_idx,$u2_idx,$h2_idx]))
    p4_gain = @lift($(cbAngAvg.checked) ? log10.(GainMatrix4Avg[2:end-1,$p1_idx,$p2_idx]) : log10.(GainMatrix4[2:end-1,$u4_idx,$h4_idx,$p1_idx,$u1_idx,$h1_idx,$p2_idx,$u2_idx,$h2_idx]))

    # limits
    y_max = @lift(max($p1_loss,$p2_loss,-40.0)+2.0)
    y_min = @lift(max($p1_loss,$p2_loss,-40.0)-12.0)
    x_min = @lift(min($p1_val,$p2_val)-2.0)
    x_max = @lift(max($p1_val,$p2_val)+2.0)

    ax1 = Axis(fig[1:2, 1:2],xlabel=L"$\log_{10}\left(p[m_\text{Ele}c]\right)$",ylabel= L"$\log_{10} \left(p\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}t} [\text{some units}]\right)$")

    onany(x_min,x_max,y_min,y_max) do _x_min, _x_max, _y_min, _y_max
        ylims!(ax1,_y_min,_y_max)
        xlims!(ax1,_x_min,_x_max)
    end

    ax1.aspect = DataAspect()

    scatterlines!(ax1,p1_val,p1_loss, label = "$name1 Loss Spectrum")
    scatterlines!(ax1,p2_val,p2_loss, label = "$name2 Loss Spectrum")   
    scatterlines!(ax1,log10.(p3_m),p3_gain, label = "$name3 Gain Spectrum")
    scatterlines!(ax1,log10.(p4_m),p4_gain, label = "$name4 Gain Spectrum")


    # lines for what the incoming state momenta are
    vlines!(ax1,p1_val,color=:black,linestyle=:dash)
    vlines!(ax1,p2_val,color=:black,linestyle=:dot)

    if name1 == "Ele" && name2 == "Pho" && name3 == "Ele" && name4 == "Pho" # Inverse Compton
        p3_val_true = @lift(log10.(replace(IC_kernel.(p2_m[$p2_idx],p1_m[$p1_idx],p2_m[$p2_idx]+sqrt(p1_m[$p1_idx]^2+1e0^2).- sqrt.(p1_m.^2 .+1e0^2)),0.0 => NaN) .* (p1_m ./ sqrt.(p1_m.^2 .+ 1e0^2)).* p1_d))
        p4_val_true = @lift(log10.(replace(IC_kernel.(p2_m[$p2_idx],p1_m[$p1_idx],p2_m),0.0 => NaN) .* p2_d))
        lines!(ax1,log10.(p1_m),p3_val_true,linestyle=:dot,linewidth=1.0, label = "ISO Inverse Compton Ele Spectrum",alpha = @lift($(cbAnalytic.checked) ? 1.0 : 0.0),color=:black)
        lines!(ax1,log10.(p2_m),p4_val_true,linestyle=:dash, linewidth=1.0, label = "ISO Inverse Compton Pho Spectrum",alpha = @lift($(cbAnalytic.checked) ? 1.0 : 0.0),color=:black)
    end
  
    Legend(subgl2[1,1:2],ax1)

    #set_close_to!(sg1.sliders[1],31)
    #set_close_to!(sg2.sliders[1],257)

    return fig

end

function InteractiveEmissionGainLossPlot(Output::Tuple)

    GLMakie.activate!(inline=false)

    (name1,name2,name3,type,m1,m2,m3,z1,z2,z3,p1_low,p1_up,p1_grid,p1_num,u1_grid,u1_num,h1_grid,h1_num,p2_low,p2_up,p2_grid,p2_num,u2_grid,u2_num,h2_grid,h2_num,p3_low,p3_up,p3_grid,p3_num,u3_grid,u3_num,h3_grid,h3_num,Ext) = Output[1]

    p1_r = DC.bounds(p1_low,p1_up,p1_num,p1_grid);
    p1_d = DC.deltaVector(p1_r);
    p1_m = DC.meanVector(p1_r);
    p2_r = DC.bounds(p2_low,p2_up,p2_num,p2_grid);
    p2_d = DC.deltaVector(p2_r);
    p2_m = DC.meanVector(p2_r);
    p3_r = DC.bounds(p3_low,p3_up,p3_num,p3_grid);
    p3_d = DC.deltaVector(p3_r);
    p3_m = DC.meanVector(p3_r);
    
    u1_r = DC.bounds(-1.0,1.0,u1_num,u1_grid);
    u2_r = DC.bounds(-1.0,1.0,u2_num,u2_grid);
    u3_r = DC.bounds(-1.0,1.0,u3_num,u3_grid);
    u1_d = DC.deltaVector(u1_r);
    u2_d = DC.deltaVector(u2_r);
    u3_d = DC.deltaVector(u3_r);
    u1_m = DC.meanVector(u1_r);
    u2_m = DC.meanVector(u2_r);
    u3_m = DC.meanVector(u3_r);
    h1_r = DC.bounds(0.0,2.0,h1_num,h1_grid);
    h2_r = DC.bounds(0.0,2.0,h2_num,h2_grid);
    h3_r = DC.bounds(0.0,2.0,h3_num,h3_grid);
    h1_d = DC.deltaVector(h1_r);
    h2_d = DC.deltaVector(h2_r);
    h3_d = DC.deltaVector(h3_r);
    h1_m = DC.meanVector(h1_r);
    h2_m = DC.meanVector(h2_r);
    h3_m = DC.meanVector(h3_r);

    fig = Figure(size = (1000,600))
    
    gl1 = GridLayout(fig[3,1:4],tellwidth = false)
    gl2 = GridLayout(fig[1:2,3:4],tellwidth = false)
    subgl1 = GridLayout(gl1[1,1],tellwidth = false)
    subgl2 = GridLayout(gl2[1,1],tellwidth = false)

    sg1 = SliderGrid(subgl1[1,1],
    (label = "p1", range = 1:length(p1_m), startvalue = 1, update_while_dragging = false),
    (label = "u1", range = 1:length(u1_m), startvalue = 1, update_while_dragging = false),
    (label = "h1", range = 1:length(h1_m), startvalue = 1, update_while_dragging = false)
    )
    sg2 = SliderGrid(subgl1[1,2],
    (label = "u2", range = 1:length(u2_m), startvalue = 1, update_while_dragging = false),
    (label = "h2", range = 1:length(h2_m), startvalue = 1, update_while_dragging = false)
    )
    sg3 = SliderGrid(subgl1[1,3],
    (label = "u3", range = 1:length(u3_m), startvalue = 1, update_while_dragging = false),
    (label = "h3", range = 1:length(h3_m), startvalue = 1, update_while_dragging = false)
    )
    sg4 = SliderGrid(subgl1[1,4],
    (label = "Ext", range = 1:length(Ext), startvalue = 1, update_while_dragging = false),
    )

    cbAngAvg = Checkbox(subgl2[2,1],checked=false)
    Label(subgl2[2,2],"Angle Averaged",halign=:left)

    p1_idx = sg1.sliders[1].value
    u1_idx = sg1.sliders[2].value
    h1_idx = sg1.sliders[3].value
    u2_idx = sg2.sliders[1].value
    h2_idx = sg2.sliders[2].value
    u3_idx = sg3.sliders[1].value
    h3_idx = sg3.sliders[2].value

    Ext_idx = sg4.sliders[1].value

    p1_val = @lift(log10.(p1_m[$p1_idx]))
    u1_val = @lift(u1_m[$u1_idx])
    h1_val = @lift(h1_m[$h1_idx])

    u3_val = @lift(u3_m[$u3_idx])
    h3_val = @lift(h3_m[$h3_idx])

    label = @lift("p1: [$(round(p1_r[$p1_idx], sigdigits=3)),$(round(p1_r[$p1_idx+1], sigdigits=3))], u1: [$(round(u1_r[$u1_idx],sigdigits=2)),$(round(u1_r[$u1_idx+1],sigdigits=2))], h1: [$(round(h1_r[$h1_idx],sigdigits=2)),$(round(h1_r[$h1_idx+1],sigdigits=2))] 
    \n 
    p2: [u2: [$(round(u2_r[$u2_idx],sigdigits=2)),$(round(u2_r[$u2_idx+1],sigdigits=2))], h2: [$(round(h2_r[$h2_idx],sigdigits=2)),$(round(h2_r[$h2_idx+1],sigdigits=2))]
    \n
    u3: [$(round(u3_r[$u3_idx],sigdigits=2)),$(round(u3_r[$u3_idx+1],sigdigits=2))], h3: [$(round(h3_r[$h3_idx],sigdigits=2)),$(round(h3_r[$h3_idx+1],sigdigits=2))]
    \n 
    Ext: [$(round(Ext[$Ext_idx],sigdigits=2))]
    ") 

    Label(subgl2[4,1:2], label)

    GainMatrix3 = Output[2]
    GainMatrix3Avg = dropdims(sum(GainMatrix3,dims=(2,3,5,6)),dims=(2,3,5,6)) / (u1_num*h1_num*u3_num*h3_num)

    # Gain Matrix
    p3_gain = @lift($(cbAngAvg.checked) ? log10.(GainMatrix3Avg[1:end,$p1_idx,$Ext_idx]) : log10.(GainMatrix3[1:end,$u3_idx,$h3_idx,$p1_idx,$u1_idx,$h1_idx,$Ext_idx]))

    # limits
    #y_max = @lift(max(maximum($p3_gain),-40.0)+2.0)
    #y_min = @lift(max(maximum($p3_gain),-40.0)-12.0)
    x_min = min(maximum(p1_r),maximum(p3_r))-1.0
    x_max = max(minimum(p1_r),minimum(p3_r))+1.0
    
    ax1 = Axis(fig[1:2, 1:2],xlabel=L"$\log_{10} \left(p[m_\text{Ele}c]\right)$",ylabel= L"$\log_{10}\left(p\frac{\mathrm{d}N}{\mathrm{d}p\mathrm{d}t} [\text{some units}]\right)$")

    onany(x_min,x_max#= ,y_min,y_max =#) do _x_min, _x_max#= , _y_min, _y_max =#
        #ylims!(ax1,_y_min,_y_max)
        xlims!(ax1,_x_min,_x_max)
    end

    ax1.aspect = DataAspect()

    #scatterlines!(ax1,p1_val,p1_loss, label = "$name1 Loss Spectrum")
    #scatterlines!(ax1,log10.(p2_m),p2_gain, label = "$name2 Gain Spectrum")  
    scatterlines!(ax1,log10.(p3_m),p3_gain, label = "$name3 Gain Spectrum")
    
    # lines for what the incoming state momenta are
    vlines!(ax1,p1_val,color=:black,linestyle=:dash)

    if type == "Sync"
        # critical frequency
        ω0 = @lift(abs((z1*1.6e-19*Ext[$Ext_idx]))/(10^$p1_val*9.11e-31))
        pc = @lift(log10(1.054e-34*$ω0/(9.11e-31*3e8^2)*(10^$p1_val)^3))
        vlines!(ax1,pc,color=:black,linestyle=(:dot,:dense))

        pstar = @lift(log10(1.054e-34*abs(z1)*1.6e-19*Ext[$Ext_idx]/((9.11e-31*3e8)^2*(sqrt((10^$p1_val)^2+m1^2)-10^$p1_val*$u1_val*$u3_val)))) 
        vlines!(ax1,pstar,color=:black,linestyle=(:dashdot,:dense))
        
    end

    return fig

end
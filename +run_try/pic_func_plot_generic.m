%% Generic plot function for oxygen - input timestep, varname,, ylabel title


function ax1 = pic_func_plot_generic(jj,varnametxt,titletxt,make_new_figure) 

    txtfile = ['./fields-' sprintf('%05d',jj) '.dat'];

    [xe,ze,ex,ey,ez,bx,by,bz,dni,dne,jix,jiy,jiz,...
    jex,jey,jez,vix,viy,viz,ti,te,nnx,nnz, ...
    dni_h,dne_h,jix_h,jiy_h,jiz_h,jex_h,jey_h,jez_h,vix_h,viy_h,viz_h,ti_h,te_h,a, ...
    wpewce,mass,pxxi_h,pyyi_h,pzzi_h,pxxe_h,pyye_h,pzze_h,vex_h,vey_h,vez_h,vex,vey,vez] ...
    = func_read_file_fields_oxygen(txtfile);


if make_new_figure == 1
    figure;
end
varname = eval(varnametxt);
ax1 = axes;
 imagesc(ax1,xe,ze,varname');
            set(gca,'YDir','normal')
            %shading interp
            cb = colorbar;
           % cb.Label.String = ylabeltxt; cb.FontSize=12;
            xlabel('X [c/\omega_{pi}]')
            ylabel('Z')
            hold on;
       ax2 = contour(xe, ze, a',20,'color',[0.3 0.3 0.3]);
            caxis([min(varname(:)) max(varname(:))]);
            figuresize(30,15,'centimeters')
            title([titletxt '       time = ' num2str(jj/wpewce(1)/mass(1)) ]);
            set(gca,'layer','top') % bring back ticks
            ylim([-15 15])

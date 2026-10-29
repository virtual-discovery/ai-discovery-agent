close all


%% GENERATE MAGNETIC STRUCTURE

spinwobj=spinw;
spinwobj.genlattice('lat_const',[5.357 5.357 13.127],'angled',[90 90 120],'spgr','R -3 m');

spinwobj.addatom('r',[1/2 0 1/2],'S',1,'label','Co2','color','b','formfact','MCo2+');

% generate canting angle for starting magnetic structure :  
% ll=5 (66 deg) ; ll=7 (72.2 deg) ; ll=10 (77 deg) ; ll=20 (83 deg) ; ll=50 (87 deg)
ll=10;
ss=[1 2 ll; 1 2 ll; 1 2 ll; -2 -1 ll; -2 -1 ll; -2 -1 ll; 1 -1 ll; 1 -1 ll; 1 -1 ll]';
ss=0.25*ss/norm(ss(:,1)');
spinwobj.genmagstr('mode','direct','k',[0 0 0],'S',ss,'nExt',[1 1 1],'unit','lu');

figMagStructure = swplot.figure;
set(gcf,'position',[10,30,300,460]);
plot(spinwobj,'range',[2 2 2],'bondMode','none');
view(2,20); zoom(2.5);

% these will be used in one of the later tests
realfhat=spinwobj.table('mag').realFhat;
[sAz, sEl, sR] = cart2sph(realfhat(:,1),realfhat(:,2),realfhat(:,3));

%return





%% GENERATE COUPLINGS

spinwobj.gencoupling('maxDistance',7);

jpars = num2cell([ -56 15 -3 ]);
[j1 j1dmi jc1] = deal(jpars{:});

spinwobj.addmatrix('label','J1','value',j1,'color','orange')
spinwobj.addmatrix('label','J1dmi','value',[0 0 j1dmi ; 0 0 0 ; -j1dmi 0 0]-0.4*[0 j1dmi 0 ; -j1dmi 0 0 ; 0 0 0],'color','orange')
spinwobj.addmatrix('label','Jc1','value',jc1,'color','purple')

spinwobj.addcoupling('mat','J1','bond',1);
spinwobj.addcoupling('mat','J1dmi','bond',1);
spinwobj.addcoupling('mat','Jc1','bond',3);

figMagBonds = swplot.figure;
set(gcf,'position',[10,600,660,520]);
plot(spinwobj,'range',[2 2 1],'magMode','none','bondMode','line','bondLineWidth0',2);

spinwobj.energy

%return





%%

% ============ OPTIMIZE MAG STRUCTURE WITH THIS HAMILTONIAN ============




%% Optimize the magnetic structure using the steep method (used in the paper)
optGroundStateSteep=true;
if optGroundStateSteep
    
    energyStart = spinwobj.energy
    
    swplot.activefigure(figMagStructure)
    plot(spinwobj,'range',[2 2 0.99],'bondMode','none');
%     pause
    optmagsteepresult = spinwobj.optmagsteep('nRun',30,'plot',true);
    
    energyOptMagStr = spinwobj.energy
    
end

%return


%% Find optimal magnetic structure through direct calculations of the energy per spin (confirms the paper result)
testStructureEnergies=false;
if testStructureEnergies
    
    spinwstructtest = spinwobj.copy
    
    
    ll=100;
    ss=[1 2 ll; 1 2 ll; 1 2 ll; -2 -1 ll; -2 -1 ll; -2 -1 ll; 1 -1 ll; 1 -1 ll; 1 -1 ll]';
    ss=0.25*ss/norm(ss(:,1)');
    spinwstructtest.genmagstr('mode','direct','k',[0 0 0],'S',ss,'nExt',[1 1 1],'unit','lu');
    energyOptMagStr = spinwstructtest.energy
    
    figMagStructure = swplot.figure;
    set(gcf,'position',[200,30,300,460]);
    plot(spinwstructtest,'range',[2 2 2],'bondMode','none');
    view(2,20); zoom(2.5);
    
    
    ll=35;
    ss=[1 2 ll; 1 2 ll; 1 2 ll; -2 -1 ll; -2 -1 ll; -2 -1 ll; 1 -1 ll; 1 -1 ll; 1 -1 ll]';
    ss=0.25*ss/norm(ss(:,1)');
    spinwstructtest.genmagstr('mode','direct','k',[0 0 0],'S',ss,'nExt',[1 1 1],'unit','lu');
    energyOptMagStr = spinwstructtest.energy
    
    figMagStructure = swplot.figure;
    set(gcf,'position',[300,30,300,460]);
    plot(spinwstructtest,'range',[2 2 2],'bondMode','none');
    view(2,20); zoom(2.5);
    
    
    ll=5;
    ss=[1 2 ll; 1 2 ll; 1 2 ll; -2 -1 ll; -2 -1 ll; -2 -1 ll; 1 -1 ll; 1 -1 ll; 1 -1 ll]';
    ss=0.25*ss/norm(ss(:,1)');
    spinwstructtest.genmagstr('mode','direct','k',[0 0 0],'S',ss,'nExt',[1 1 1],'unit','lu');
    energyOptMagStr = spinwstructtest.energy
    
    figMagStructure = swplot.figure;
    set(gcf,'position',[400,30,300,460]);
    plot(spinwstructtest,'range',[2 2 2],'bondMode','none');
    view(2,20); zoom(2.5);
    
    
    
    ll=0; %5->66deg %10000;
    ss=[1 2 ll; 1 2 ll; 1 2 ll; -2 -1 ll; -2 -1 ll; -2 -1 ll; 1 -1 ll; 1 -1 ll; 1 -1 ll]';
    ss=0.25*ss/norm(ss(:,1)');
    spinwstructtest.genmagstr('mode','direct','k',[0 0 0],'S',ss,'nExt',[1 1 1],'unit','lu');
    energyOptMagStr = spinwstructtest.energy
    
    figMagStructure = swplot.figure;
    set(gcf,'position',[500,30,300,460]);
    plot(spinwstructtest,'range',[2 2 2],'bondMode','none');
    view(2,20); zoom(2.5);

end

% return


%% Find optimal magnetic structure through the generic 3D relaxation (not perfect, but shows that the previous two are robust)
optGroundStateConstrained=false;
if optGroundStateConstrained
    

    % energy of the known helical configuration w.r.t. these coupling parameters
    energyStart = spinwobj.energy;
    
    % starting configuration (angles from definition of the spin structure above)
    sAngles=[sEl ; sAz];
    
    optRes = spinwobj.optmagstr('func',@gm_spherical3dd,'nRun',3,'x0',sAngles(:)');
    
    energyOptMagStr = optRes.e;
    
    % plot new structure
    hFigure = swplot.figure;
    % plot(spinwobj,'range',[1 1 1],'bondMode','line','bondLineWidth0',2)
    plot(spinwobj,'range',[2 2 1],'bondMode','none');

end

%return



%%

% ============ MAKE DISPERSION PLOTS ============


%% Plot the spin wave dispersion along various directions
makeDispersionPlot=false;
if makeDispersionPlot
    
    Qlist = {[0 0 0] [0 0 3] [3 0 3] [3 3 3] [0 0 3] 100};
    Qlab  = {'(000)' '(003)' '(303)' '(333)' '(003)'};
    
    coSpec = spinwobj.spinwave(Qlist,'hermit',false);
    coSpec = sw_neutron(coSpec);
    coSpec = sw_egrid(coSpec,'component','Sperp','imagChk',false);
    figure
    set(gcf,'position',[400,600,750,430])
    sw_plotspec(coSpec,'mode',3,'axLim',[0 5],'dE',10,'dashed',true,'qlabel',Qlab)
    colorbar off
    
    
end


%% Plot the spin wave dispersion showing the low-energy gap
makeGapPlot=false;
if makeGapPlot
    
    Qlist = {[0 0 0] [0 0 3] [0 0 6] 30};
    Qlab  = {'(000)' '(003)' '(006)' };
    
    gapSpec = spinwobj.spinwave(Qlist,'hermit',false);
    
    gapSpec = sw_neutron(gapSpec);
    gapSpec = sw_egrid(gapSpec,'component','Sperp','imagChk',false);
    figure
    set(gcf,'position',[200,600,750,430])
    sw_plotspec(gapSpec,'mode',3,'axLim',[0 25],'dE',1,'dashed',true,'qlabel',Qlab)
    ylim([0 25])
    colorbar off
    
    annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'string', jpars)
    
    
end






%%

% ============ MAKE 2D SLICE PLOTS ============



%% Make constant energy slices using Horace, with IN8 data points overlayed
%horace_on()
makeConstEnergyPlots=false;
if makeConstEnergyPlots
    
    dq=0.025; % choose the Q resolution here for all cuts (faster plotting when resolution is coarse)
    
    Ebin   = [14,3,20];
    fwhm0  = 3;
    d3dobj = d3d(spinwobj.abc,[1 0 0 0],[-3,dq,3],[0 0 1 0],[-9,dq,0],[0 0 0 1],Ebin);
    d3dobj = disp2sqw_eval(d3dobj,@spinwobj.horace,{'component','Sperp','formfact',true,'imagChk',false},fwhm0);
    d3dobj = smooth(d3dobj);
    plot(cut(d3dobj,[],[],[Ebin(1) Ebin(3)]));
    set(gcf,'position',[680,700,370,460])
    colorslider('delete'); set(gca,'YDir','normal'); axis square; axis tight; colorbar; colormap(jet);
    hold on
    plot([0.1000,0.1040],-[3.0000,6.0000],'r*')
    hold off
    keep_figure

    Ebin   = [29,3,35];
    fwhm0  = 3;
    d3dobj = d3d(spinwobj.abc,[1 0 0 0],[-3,dq,3],[0 0 1 0],[-9,dq,0],[0 0 0 1],Ebin);
    d3dobj = disp2sqw_eval(d3dobj,@spinwobj.horace,{'component','Sperp','formfact',true,'imagChk',false},fwhm0);
    d3dobj = smooth(d3dobj);
    plot(cut(d3dobj,[],[],[Ebin(1) Ebin(3)]));
    set(gcf,'position',[1050,700,370,460])
    colorslider('delete'); set(gca,'YDir','normal'); axis square; axis tight; colorbar; colormap(jet); %caxis([0 3])
    hold on
    plot([0.1770,0.1750,0.1880,0.1850,0.1990,0.1540,0.1390,0.1760,0.1570,0.0891],-[6.0000,6.0000,6.2000,6.4000,6.6000,6.8000,7.0000,6.4000,6.8000,7.2000],'r*')
    hold off
    keep_figure
    
    Ebin   = [47,5,57];
    fwhm0  = 3;
    d3dobj = d3d(spinwobj.abc,[1 0 0 0],[-3,dq,3],[0 0 1 0],[-9,dq,0],[0 0 0 1],Ebin);
    d3dobj = disp2sqw_eval(d3dobj,@spinwobj.horace,{'component','Sperp','formfact',true,'imagChk',false},fwhm0);
    d3dobj = smooth(d3dobj);
    plot(cut(d3dobj,[],[],[Ebin(1) Ebin(3)]));
    set(gcf,'position',[1420,700,370,460])
    colorslider('delete'); set(gca,'YDir','normal'); axis square; axis tight; colorbar; colormap(jet); %caxis([0 3])
    hold on
    plot([0.2370,0.2530,0.2130,0.2290,0.2190],-[6.0000,6.2000,6.4000,5.8000,5.6000],'r*')
    hold off
    keep_figure


end




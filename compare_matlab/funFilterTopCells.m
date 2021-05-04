function CBfilteredBool=funFilterTopCells(M1,nTopCells)

%% filtered CB: top 5000 cells
umiPerCB=full(sum(M1,1));
umiPerCBsorted=sort(umiPerCB,'descend');
CBfilteredBool=umiPerCB>=umiPerCBsorted(nTopCells);


function myhandle =pattern_extrand(strip)
cnt1=0;
cnt2=strip/2;
%pitch=1;
for i=1:1
    if i==1
   myhandle= plot([cnt1 (cnt1+(strip/2))], [0 0], '-c', 'LineWidth',2);
    hold on
    else
      myhandle= plot([cnt1 (cnt1+strip)], [0 0], '-c', 'LineWidth',2);
      hold on
    end
 myhandle =plot([cnt2 (cnt2+strip)], [0 0], '-k', 'LineWidth',2);
hold on
    
    if i==1
        cnt1=cnt1+1.5*strip;
    else
        cnt1=cnt1+2*strip;
    end
    cnt2=cnt2+2*strip;
end
cnt1=0;
cnt2=-strip/2;
for i=1:1
    if i==1
     myhandle=plot([cnt1 (cnt1-(strip/2))], [0 0], '-c', 'LineWidth',2);
    hold on
    else
      myhandle= plot([cnt1 (cnt1-strip)], [0 0], '-c', 'LineWidth',2);
      hold on
    end
plot([cnt2 (cnt2-strip)], [0 0], '-k', 'LineWidth',2);
hold on
    
    if i==1
        cnt1=cnt1-1.5*strip;
    else
        cnt1=cnt1-2*strip;
    end
    cnt2=cnt2-2*strip;
end
end
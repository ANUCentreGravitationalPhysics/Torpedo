function[ufc,lb,ab,hb,st] = opt(m,mt,smax,lmax,wmax)
%calculates the optimum blade parameters
% changed 15/4/03 to allow finer steps in thickness

g = 9.81; %gravitational acceleration
%alpha = 1.38; %shape factor
%alpha = 1.53; %shape factor alternate value
alpha = 1.35; %average for larger blades, from M Plissi Mar04
ye    = 186e9; %Young's elastic modulus for Marval steel blades
%ye=176e9;% alternative value
mtb = mt;%total load on one blade
mb  = m;%uncoupled mass on one blade
lb = lmax; %always best 
for index = 1:60
   h = 0.001 + index*0.0001; 
   for jndex = 1:18
      a = 0.03+(wmax-0.03)/18*jndex;  
      s = 6.*mtb.*g.*lb./(a.*h.^2);%max stress at equilibrium
      d = 4.*mtb.*g.*lb.^2.*alpha./(ye.*a.*h.^3); %deflection unloaded /l
      if s < smax & d < pi/2
         f(index,jndex) = sqrt(ye.*a.*h.^3./(4.*mb.*lb.^3.*alpha))./2./pi; %mode freq.
         stress(index,jndex) = s;
      else 
         f(index,jndex) = NaN;%invalidate 
         stress(index,jndex) = s;
      end
   end
end
[fmin1,index1] = min(f);
[fmin2,index2] = min(fmin1);
ufc = fmin2;
aindex = index2;
ab     = 0.03+(wmax-0.03)/18*aindex;  
hindex = index1(aindex);
hb = 0.001 + hindex*0.0001; 
st = stress(hindex,aindex);

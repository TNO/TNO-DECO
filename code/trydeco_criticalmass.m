function usemasses=trydeco_criticalmass(data,residuelimit)
% 
[evectorm,evaluem] =eig(data'*data); % eigen decomposition based on mass
mindex=deco_criticalmassindex(evectorm);
last=length(mindex);
for i=last:-1:1
    rv=deco_rvcoeff(data,data(:,mindex(i:last)));
    if(rv==0)
        msgbox('No significant information in block','warn');
        usemasses=[];
        break;
    end
    if(rv>=1-residuelimit) % cutoff has to be statistically defined; to be implemented
        usemasses=mindex(i:last);
        break;
    end
end

function y=deco_criticalmassindex(v)
[coeff,index]=sort(abs(v),'descend');
endval=size(index,2);
y=[];
for j=endval:-1:1
    i=1;
    while (i<=size(index,1))
        ty=index(i,j);
        if(find(y==ty))
            i=i+1;
        else
            y(j)=ty;
            break;
        end
    end
end


function y=deco_rvcoeff(d1,d2)
% mcd1=deco_meancenter(d1);
% mcd2=deco_meancenter(d2);

%s=mcd1*mcd1';
%t=mcd2*mcd2';
s=d1*d1';
t=d2*d2';
y=double(trace(s*t))/sqrt(double(trace(s*s))*double(trace(t*t)));


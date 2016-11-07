function Schedule=FireSchedule(mean,shape,years,delay)
% Creates a schedule of fires based on a Weibull distrobibution with
% parameters mean and shape
    if nargin<4
        delay=0;
    end
scale=mean/(gamma(1+(1/shape)));
current=0;
i=1;
IntLen=[];
IntLen(round(years/scale))=0;
    % Need a separate case so that IntLen(i)>0.5
    while(current<0.5)
        U=rand(1);
        IntLen(i)=scale*(-log(U))^(1/shape)+delay;
        current=IntLen(i);
    end
    i=i+1; 
    while current<years
        U=rand(1);
        IntLen(i)=scale*(-log(U))^(1/shape)+delay;
        current=current+IntLen(i);
        i=i+1;
    end
    Total=0;
    Schedule=zeros(1,years);
    for j=1:length(IntLen)
        Total=Total+IntLen(j);
        Schedule(round(Total))=1;
    end
    Schedule=[1 Schedule]; 
    Schedule=Schedule(1:years);
end


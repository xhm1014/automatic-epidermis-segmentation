function PS=XFilter(PI,ws,method)
if ~exist('method','var')
    method='average';
end

if strcmp(method,'average')
    w=floor(ws/2);
    X=padarray(PI(:,1),[w 0],PI(1,1),'pre');
    X=padarray(X,[w 0], PI(end,1),'post');
    
    Y=padarray(PI(:,2),[w 0],PI(1,2),'pre');
    Y=padarray(Y,[w 0],PI(end,2),'post');
    F=ones(1,ws)/ws;
    PS(:,1)=conv(X,F,'valid');
    PS(:,2)=conv(Y,F,'valid');
end
end
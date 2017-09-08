% cloud scheme wrapper function.
THETA_A=200:5:400;
P_A=20000:1000:1.1e5;
QT_A=linspace(1e-4,40e-3,30);
QC_A=zeros(length(theta),length(P),length(qt));
QV_A=zeros(length(theta),length(P),length(qt));
TC_A=zeros(length(theta),length(P),length(qt));

for i=1:length(THETA_A)
    for j=1:length(P_A)
        for k=1:length(QT_A)
            [QC_A(i,j,k),QV_A(i,j,k),TC_A(i,j,k)]= ...
                cloud_scheme(THETA_A(i),P_A(j),QT_A(k));
        end
    end
end
save cloud_scheme_data THETA_A P_A QT_A QC_A QV_A TC_A;

function beta=solvebeta(N,d,epsilon,wsd,rr)
V = 1:N;
C = nchoosek(V,N-d);
beta=0;
for i=1:size(C,1)
    beta_sub=1;
    for j=1:N-d
        %beta_sub=beta_sub*(1-epsilon);
        if (epsilon-wsd(C(i,j))/rr(C(i,j)))<0
            beta_sub=beta_sub;
        else
            beta_sub=beta_sub*(1-epsilon+wsd(C(i,j))/rr(C(i,j)));
        end
    end
    beta=beta+beta_sub;
end
end

function [ eta] = stepsize_bisection(eval_cost,eval_gradient,parameters )
grad=eval_gradient(parameters);
p=-grad;
a=0;
b=2;
c1=0.2;
c2=0.3;
goahead=1;

while(eval_cost(parameters+b*p)<=eval_cost(parameters)+c1*b*grad'*p)
        b=2*b;
end

while(goahead==1)
        eta=(a+b)/2;
        q=eval_gradient(parameters+eta*p);
        if(eval_cost(parameters+eta*p)>eval_cost(parameters)+c1*eta*grad'*p)
                b=eta;
        elseif(q'*p<c2*grad'*p)
                a=eta;
        else
                goahead=0;
        end
end


end


function [thissim,bestsimm]=calcsimm(X,test);

if exist('test','var'),
    % test on size of X and test
    [m,n]=size(X);
    [o,p]=size(test);
    X=[X; test];
    thissim=calcsimm(X);
    size(thissim);
    temp=thissim-eye(m+o);
    for i=1:o,
        [bestsimm.value(i),bestsimm.id(i)]=max(temp(:,m+i));
    end
    thissim=thissim(1:m,1:m);
else
    [m,n]=size(X);
    if m==1,
        thissim=1;
    else
        these=nchoosek([1:m],2); % find all possible combinations of rows and cols
        thissim=eye(m); % initialize diagonal of matrix 
        for i=1:size(these,1),
            % normalized length of (i,1) element
            A_lambda_ref  = X(these(i,1),:)*sparse(diag((1./sqrt(sum(X(these(i,1),:).^2)))));
            % normalized length of (i,2) element
            A_lambda_test = X(these(i,2),:)*sparse(diag((1./sqrt(sum(X(these(i,2),:).^2)))));
            thissim(these(i,1),these(i,2))=sum(A_lambda_ref*A_lambda_test');
            thissim(these(i,2),these(i,1))=thissim(these(i,1),these(i,2));
        end
    end
end

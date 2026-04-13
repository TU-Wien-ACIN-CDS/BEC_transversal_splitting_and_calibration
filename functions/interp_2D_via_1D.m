function S_x = interp_2D_via_1D(X,S,Xv)


N_x = numel(Xv);
[N_y,~] = size(S);

S_x = zeros(N_y,N_x);

% thresh = 1e-9;
% 
% for i = 1:N_x
%     xv = Xv(i);
%     if any(abs(xv-X)<thresh)
%         S_x(:,i) = S(:,abs(xv-X)<thresh);
%     elseif xv > X(end)
%         S_x(:,i) = S(:,end);
%     elseif xv < X(1)
%         S_x(:,i) = S(:,1);
%     else
%         indx_pre    = find(X<xv,1,"last");
%         indx_post   = find(X>xv,1,"first");
%         X_pre       = X(indx_pre);
%         X_post      = X(indx_post);
%         S_x(:,i)    = 1/(X_post - X_pre)*((xv - X_pre)*S(:,indx_pre)+(X_post - xv)*S(:,indx_post));
%     end
% end

for i = 1:N_y
    S_x(i,:) = interp1(X,S(i,:),Xv,"cubic");
end

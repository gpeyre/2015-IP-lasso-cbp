function install_cvx(cvx_path)

if nargin==0
    cvx_path = '/Users/gpeyre/Dropbox/work/sparsity-superresolution/codes/cvx';
end
if not(exist('cvx_startup'))
    p = pwd;
    cd(cvx_path);
    cvx_startup;
    cd(p);
end
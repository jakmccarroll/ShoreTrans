%% figpos -> set position of figure
% function figpos(ox,oy,sx,sy) 
    % set(gcf,'position',[sz(3)*ox sz(4)*oy sz(3)*sx sz(4)*sy]);
    % INPUTS
        % ox, oy -> fig x,y origin as a ratio of screensize
        % sx,sy -> fig x,y size as a ratio of screensize

function figpos(varargin)

if nargin==0
    ox=.3;
    oy=.2;
    sx=.5;
    sy=.6;
else
    ox=varargin{1};
    oy=varargin{2};
    sx=varargin{3};
    sy=varargin{4};
end
    
sz=get(0,'screensize');
set(gcf,'position',[sz(3)*ox sz(4)*oy sz(3)*sx sz(4)*sy]);

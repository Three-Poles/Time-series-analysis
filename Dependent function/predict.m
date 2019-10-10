function [yhat,xi,thpred]=predict(varargin)
%PREDICT Computes the k-step ahead prediction.
%   YP = PREDICT(MODEL,DATA,K)
%
%   DATA: The output - input data as an IDDATA object, for which the
%   prediction is computed.
%
%   MODEL: The model as any IDMODEL object, IDPOLY, IDSS, IDARX or IDGREY.
%
%   K: The prediction horizon. Old outputs up to time t-K are used to
%       predict the output at time t. All relevant inputs are used.
%       K = Inf gives a pure simulation of the system.(Default K=1).
%   YP: The resulting predicted output as an IDDATA object. If DATA
%       contains multiple experiments, so will YP.
%
%   YP = PREDICT(MODEL,DATA,K,INIT)  or
%   YP = PREDICT(MODEL,DATA,K,'InitialState',INIT) allows the choice of
%   initial state vector:
%   INIT: The initialization strategy: one of
%      'e': Estimate initial state so that the norm of the
%           prediction errors is minimized.
%           This state is returned as X0e, see below. For multiexperiment
%           DATA, X0e is a matrix whose columns contain the initial states
%           for each experiment.
%      'z': Take the initial state as zero
%      'm': Use the model's internal initial state.
%      'd': Same as 'e', but if Model.InputDelay is non-zero, these delays
%            are first converted to explicit model delays, so that the are
%            contained in X0 for the calculation of YP.
%       X0: a column vector of appropriate length to be used as initial value.
%           For multiexperiment DATA, X0 may be a matrix whose columns give
%           different initial states for each experiment.
%   If INIT is not specified, Model.InitialState is used, so that
%      'Estimate', 'Backcast' and 'Auto' gives an estimated initial state,
%      while 'Zero' gives 'z' and 'Fixed' gives 'm'. If Model is an IDARX
%      model the default initial state is 'z'.
%
%   With [YP,X0e,MPRED] = PREDICT(MODEL,DATA,K) the initial state(s) and the
%   predictor MPRED are returned.  Note that if MODEL is continuous time, X0e
%   is returned as the states of this model. These may differ from the initial
%   states of the discrete time model that has the sampling interval(s) of
%   DATA. Also when INIT = 'd' only the states of MODEL are returned in
%   X0e. (To obtain the full X0, apply first INPD2NK to a discrete time
%   version of the model.)
%   MPRED is a cell array of IDPOLY objects,
%   such that MPRED{ky} is the predictor for the ky:th output. The matching
%   the channels of MPRED with data follows from its InputNames.
%   See also COMPARE and IDMODEL/SIM.

%   L. Ljung 10-1-89,9-9-94
%   Copyright 1986-2007 The MathWorks, Inc.
%   $Revision: 1.1.8.7 $  $Date: 2007/12/14 14:44:34 $


% First find if there is any pair 'InitialState', init in the argument
% list:
nr = find(strncmpi(varargin,'in',2));
init =[]; xi = [];
if ~isempty(nr)
    if length(varargin)<nr+1
        error('ident:analysis:predictCheck1',...
            '''InitialState'' must be followed by a character or a vector/matrix.')
    end
    init = varargin{nr+1};
    if ~isa(init,'char') && ~isa(init,'double')
        error('ident:analysis:predictCheck2',...
            'The initial state must either be a character or a vector/matrix.')
    end
    varargin(nr:nr+1)=[];%
end
if length(varargin)<2
    disp('Usage: YP = PREDICT(MODEL,DATA)')
    disp('       YP = PREDICT(MODEL,DATA,M)')
    yhat = []; xi=[];
    return
end
data = varargin{1};
theta = varargin{2};
if length(varargin)>2
    m = varargin{3};
    errm = 0;
    if ~isnumeric(m)
        errm = 1;
    end
    if round(m)~=m || m<1
        errm = 1;
    end
    if errm
        error('ident:analysis:predictCheck3',...
            'The prediction horizon must be a positive integer.')
    end
else
    m = 1;
end
if length(varargin)>3 && isempty(init)
    init = varargin{4};
end
thpred = [];
yhat = [];
if isa(data,'idmodel') % Forgive order
    data1 = theta;
    theta = data;
    data = data1;
end
nxorig = length(ssdata(theta));
[nym,num] = size(theta);
if isa(data,'frd') || isa(data,'idfrd')
    error('ident:analysis:predictCheck4',...
        '"predict" command is not applicable to frequency response data (IDFRD).')
end
if isa(data,'iddata')
    [Ncap,ny,nu,Ne] = size(data);
    if ny~=nym || nu~=num
        error('ident:estimation:modelDataDimMismatch',...
            'The input/output size of the model does not match with that of data.')
    end
end

if isempty(m), m=1; end
if ~isinf(m),
    if m<1 || m~=floor(m)
        error('ident:analysis:predictCheck3',...
            'Prediction horizon must be a positive integer.')
    end
end
% First deal with the special case of continuous time data and continuous
% model:
if isa(data,'iddata') && any(cell2mat(pvget(data,'Ts'))==0)
    if ~isinf(m)
        error('ident:analysis:predictCTData',...
            'The "predict" command cannot be used for continuous time data.')
    end
    yhat = sim(data,theta);
    if nargout == 0
        utidplot(theta,yhat,'Simulated')
        clear yhat
    end
    return
end

if isempty(init)
    if isa(theta,'idarx')
        init='z';
    else
        init=pvget(theta,'InitialState');
    end
    init=lower(init(1));
    if init=='d'
        inpd= pvget(theta,'InputDelay');
        if norm(inpd)==0
            init = 'e';
        end
    end
    if init=='a' || init=='b'
        init = 'e';
    elseif init=='f'
        init = 'm';
    end
elseif ischar(init)
    init=lower(init(1));
    if ~any(strcmp(init,{'m','z','e','d'}))
        error('ident:analysis:predictCheck5',...
            'The fourth input, "INIT", to "predict command" must be one of ''E(stimate)'', ''Z(ero)'', ''M(odel)'', ''D(elayconvert)'', or a vector.')
    end
    if init=='m'&&isa(theta,'idpoly')
        warning('ident:analysis:predictCheck6','For an IDPOLY Model, ''INIT'' = ''m'' implies ''Zero'' state values.')
    end
end

%% Now init is either a vector or the values 'e', 'm', 'z' 'd'
if isa(data,'iddata');
    if strcmp(pvget(data,'Domain'),'Frequency')
        if m<inf
            error('ident:analysis:predictCheck7',...
                'Prediction with finite horizon is not applicable to frequency domain IDDATA.')
        end
    end
end
if isa(data,'iddata')
    [uni,Tsd,inters] = dunique(data);
    iddflag = 1;
else
    iddflag = 0;
    uni = 1;
    Tsd = 1;
    inters = 'z';
end
%% Do away with continuous time right away:
Ts = pvget(theta,'Ts');
if Ts == 0
    theta = pvset(theta,'CovarianceMatrix',[]);% Less calculations

    ms = idss(theta);
    nc = size(ms,'nx');
    if ischar(init)
        if init=='m';
            init=pvget(ms,'X0');
        end
        if init=='z'
            init=zeros(nc,1);
        end
    end
    if ~isa(data,'iddata')
        error('ident:analysis:predictCheck8',...
            ['For a continuous time model, data must be given as an IDDATA object.',...
            '\nUse DATA = iddata(y,u,Ts).'])

    end
    Nex = size(data,'Ne');

    u=pvget(data,'InputData');
    u1 = [];
    for kexp = 1:length(u);
        u1=[u1,u{kexp}(1,:).'];
    end
    nu = size(data,'nu');

    if uni
        nc = size(ms,'Nx');
        [md,Gll] = c2d(ms,Tsd,inters);
        if init=='d', md = inpd2nk(md); init='e';end

        % First work on init
        if isa(init,'double')
            [nxi,mexp]=size(init);
            if nxi~=nc
                error('ident:analysis:predictCheck9','The initial states should have %d row(s).',nc)
            end
            if mexp==1 && Nex>1
                init = init*ones(1,Nex);
                mexp = Nex;
            end
            if mexp~=Nex
                error('ident:analysis:predictCheck10',...
                    'When the initial state is specified, the number of columns must either be 1 or equal to the number of experiments.');
            end
            if ~isempty(Gll)
                init = Gll*[init;u1];
            end

        end
        if init == 'e'
            inite = x0iniest(md,data);
        else
            inite = init;
        end
        yhat = predict(md,data,m,inite); %Check XIC!
        if ~isempty(Gll)
            xi = inite(1:nc,:)-Gll(1:nc,nc+1:nc+nu)*u1;
        else
            xi = inite(1:nc,:);
        end

    else % Different sampling intervals
        Tsd = pvget(data,'Ts');
        ints = pvget(data,'InterSample');
        xic=zeros(nc,Nex);
        for kexp = 1:Nex
            [md,Gll] = c2d(theta,Tsd{kexp},ints{1,kexp});
            if init == 'd', md = inpd2nk(md); init = 'e';end
            if isa(init,'double')
                if ~isempty(Gll)
                    if size(init,2)==1
                        initk = Gll*[init;u1(:,kexp)];
                    else
                        if size(init,2)~=Nex
                            error('ident:analysis:predictCheck10',...
                                'When the initial state is specified, the number of columns must either be 1 or equal to the number of experiments.');
                        end
                        initk = Gll*[init(:,kexp);u1(:,kexp)];
                    end
                end

            else
                initk = init;
            end
            if init == 'e'
                initk = x0iniest(md,getexp(data,kexp));
            end
            if ~isempty(Gll)
                xi(:,kexp) = initk(1:nc)-Gll(1:nc,nc+1:nc+nu)*u1(:,kexp);
            else
                xi(:,kexp) = initk(1:nc);
            end
            yhatk = predict(md,getexp(data,kexp),m,initk);
            if isempty(yhat)
                yhat = yhatk;
            else
                yhat=merge(yhat,yhatk);
            end
        end
    end
    if nargout == 0
        utidplot(md,yhat,'Predicted')
        clear yhat
    end
    return
else %Ts>0
    if init=='d'
        if any(pvget(theta,'InputDelay')>0)
            theta = inpd2nk(theta);
        end
        init = 'e';
    end
end

if isinf(m)
    if ischar(init) && init=='e'
        X0 = x0iniest(theta,data);

        init = X0;  % Note the difference with COMPARE, which fits
        % X0 with K = 0;
    end
    if ~isa(data,'iddata')
        data = data(:,nym+1:end);
    end
    yhat = sim(theta,data,init);
    thpred = theta;
    if nargout == 0
        utidplot(theta,yhat,'Predicted')
        clear yhat
    elseif nargout>1
        xi = init;
    end
    return
end
Inpd = pvget(theta,'InputDelay');
nu = size(theta,'nu');

if ~uni
    Tsd = pvget(data,'Ts'); Tsd = Tsd{1};
    warning('ident:estimation:nonUniqueDataTs',...
        ['The data set contains experiments with different sampling intervals.\n',...
        'Sampling interval from the first experiment (=%g) will be used.'],Tsd)
end

if iddflag && any(abs(Ts-Tsd)>1e4*eps)
    warning('ident:analysis:predictCheck11',...
        'The data and model sampling intervals are different.')
end

Inpd = pvget(theta,'InputDelay');
if init=='d'
    init = 'e';

    theta = inpd2nk(theta);
    Inpd = zeros(size(Inpd));
end
if isa(data,'iddata')
    if isnan(data)
        error('ident:analysis:missingData',...
            'IDDATA object contains missing data.\nUse "misdata" command to fill in missing (NaN) samples before using the data object.');
    end
    data = nkshift(data,Inpd,'append');
    theta = pvset(theta,'InputDelay',zeros(nu,1));
    [ze,Ne,ny,nu,Ts,Name,Ncaps,errflag] = idprep(data,0,'dummy');
    %if ~isempty(errflag.message), error(errflag), end
    if ~isempty(Name), data.Name = Name; end
    iddatflag = 1;
else
    if norm(Inpd)>eps

        if iscell(data)  
            data = data{1};
        end
        [Ncap,nudum] = size(data);
        nk1 = Inpd;
        ny = nudum - length(nk1);

        Ncc = min([Ncap,Ncap+min(nk1)]);
        for ku = 1:length(nk1)
            u1 = data(max([nk1(:);0])-nk1(ku)+1:Ncc-nk1(ku),ny+ku);
            newsamp = Ncap-length(u1);
            if nk1(ku)>0
                u1= [zeros(newsamp,1);u1];
            else
                u1 = [u1;zeros(newsamp,1)];
            end
            data(:,ny+ku) = u1;
        end

    end
    iddatflag= 0;
    if ~iscell(data), ze ={data};else ze = data;end
    Ne = length(ze);
    nz = size(ze{1},2);
    for kexp = 1:Ne
        Ncaps(kexp) = size(ze{kexp},1);
    end
    Ts =[];

end
if m>=min(Ncaps)
    error('ident:analysis:predictCheck12',...
        'Prediction horizon must be less than the number of data samples.')
end

[A,B,C,D,K,X0]=ssdata(theta);
[nyt,nx]=size(C);nut=size(B,2);
if iddatflag
    if ny~=nyt || nu~=nut
        error('ident:estimation:modelDataDimMismatch',...
            'The input/output size of the model does not match with that of data.')
    end
else
    if nz~=nyt+nut
        error('ident:estimation:modeldoubleDataDimMismatch',...
            'Data matrix size inconsistent with number of model input/output.')
    end
    ny = nyt; nu = nut;
end
if strcmp(init,'e')
    X0 = x0iniest(theta,ze);

elseif strcmp(init,'z')
    X0 = zeros(size(X0));
end
if ~ischar(init)
    X0 = init;
end
if nargout>1
    xi = X0;
end
[xnr,xnc]= size(X0);
if xnc~=1 && xnc~=Ne
    error('ident:analysis:predictCheck10',...
        'When the initial state is specified, the number of columns must either be 1 or equal to the number of experiments.');

end
if xnr~=nx
    error('ident:analysis:predictCheck13','Initial state matrix must have %d rows.',nx)
end

if xnc==1 && Ne>1
    X0= X0*ones(1,Ne);
end
for kexp = 1:Ne
    z = ze{kexp};
    u = z(:,1+ny:end);
    Ncap = Ncaps(kexp);

    if m==inf,
        yhat=sim(theta,u,X0(:,kexp));
    else
        x=ltitr(A-K*C,[K B-K*D], z, X0(:,kexp));
        if m==1,
            yhat=(C*x.').';
            if ~isempty(D),yhat=yhat + (D*u.').';end
        else
            F=D;Mm=eye(length(A));
            for km=1:m-1
                F=[F C*Mm*B];
                Mm=A*Mm;
            end
            yhat=zeros(Ncap,ny);%corr 911111
            for ky=1:ny
                for ku=1:nu
                    yhat(:,ky)=yhat(:,ky)+filter(F(ky,ku:nu:m*nu),1,u(:,ku));
                end
            end
            if isempty(yhat),yhat=zeros(Ncap,ny);end
            yhat(m:Ncap,:)=yhat(m:Ncap,:)+(C*Mm*x(1:Ncap-m+1,:).').';
            if nu>0
                x=ltitr(A,B,u(1:m,:),X0(:,kexp));
                yhat(1:m,:)=(C*x.').';
            end
            if ~isempty(D),yhat(1:m,:)=yhat(1:m,:) + (D*u(1:m,:).').';end
        end
    end
    yhatc{kexp} = yhat;
end
if iddatflag
    yhat = data;
    yhat = pvset(yhat,'OutputData',yhatc,'InputData',[]);
else
    yhat = yhatc;
end
if nargout >2
    thpred = polypred(theta,m);
end
if nargout > 1
    try
        xi = xi(1:nxorig,:);
    end
end
if (nargout == 0)
    % Plot y and yhat.
    utidplot(theta,yhat,'Predicted');
    clear yhat x0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xi = x0iniest(m,data)

if isa(data,'iddata') && strcmpi(pvget(data,'Domain'),'frequency')
    [dum,xi] = pe_f(m,data);
    return
end

Inpd = pvget(m,'InputDelay');
if isa(data,'iddata')
    data = nkshift(data,Inpd,'append');
    [ze,Ne,ny,nu,Ts,Name,Ncaps,errflag] = idprep(data,0,'dummy');
    %if ~isempty(errflag.message), error(errflag), end
    if ~isempty(Name)
        data.Name = Name;
    end
else
    if iscell(data)
        ze = data;
    else
        ze = {data};
    end
end

alg=pvget(m,'Algorithm');
maxsize=alg.MaxSize;

[A,B,C,D,K,X0]=ssdata(m);

nx=length(A);
if ischar(maxsize)
    maxsize = idmsize(length(ze{1}),nx);
end

AKC=A-K*C;
[ny,nx]=size(C);
nu=size(B,2);
el=zeros(0,ny);
xic = zeros(nx,0);
Ne = length(ze);
for kexp = 1:Ne
    z = ze{kexp};
    [Ncap,nz]=size(z);
    if nu+ny~=nz,
        error('ident:analysis:predictCheck14',...
            ['Incorrect number of data columns specified.',...
            '\nIt should be equal to the sum of the number of inputs and the number of outputs of the model.'])
    end

    nz=ny+nu; [Ncap,dum] = size(z); n = nx;
    rowmax = nx+nz; X0 = zeros(nx,1);
    M=floor(maxsize/rowmax);
    if ny>1 || M<Ncap
        R=zeros(n,n);Fcap=zeros(n,1);R1=[];
        for kc=1:M:Ncap
            jj=(kc:min(Ncap,kc-1+M));
            if jj(length(jj))<Ncap
                jjz = [jj,jj(length(jj))+1];
            else
                jjz=jj;
            end
            psitemp=zeros(length(jj),ny);
            psi=zeros(ny*length(jj),n);
            x=ltitr(AKC,[K B-K*D],z(jjz,:),X0);
            yh=(C*x(1:length(jj),:).').';
            nm=pvget(m,'NoiseVariance');
            if isempty(nm) || norm(nm)==0 % To handle models without noise
                nm = eye(ny);
            end
            sqrlam=pinv(sqrtm(nm));
            if ~isempty(D),yh=yh+(D*z(jj,ny+1:ny+nu).').';end

            e=(z(jj,1:ny)-yh)*sqrlam;
            [nxr,nxc]=size(x);X0=x(nxr,:).';
            evec=e(:);
            kl=1;
            for kx=1:nx
                if kc==1
                    x0dum=zeros(nx,1);x0dum(kx,1)=1;
                else
                    x0dum=X00(:,kl);
                end
                psix=ltitr(AKC,zeros(nx,1),zeros(length(jjz),1),x0dum);
                [rp,cp]=size(psix);
                X00(:,kl)=psix(rp,:).';
                psitemp=(C*psix(1:length(jj),:).').'*sqrlam;
                psi(:,kl)=psitemp(:);kl=kl+1;
            end
            if ~isempty(R1)
                if size(R1,1)<n+1
                    error('ident:estimation:X0EstTooSmallMaxSize',...
                        'The value of "MaxSize" algorithm property is too small to estimate the initial states.')
                end
                R1=R1(1:n+1,:);
            end
            H1 = [R1;[psi,evec] ];R1 = triu(qr(H1));

        end
        try
            xi(:,kexp) = pinv(R1(1:n,1:n))*R1(1:n,n+1);
        catch
            warning('ident:estimation:X0EstFailed',...
                'Failed to estimate initial conditions. Check model stability. \nInitial state has been set to zero.')
            xi(:,kexp) = zeros(n,1);
        end
    else
        %% First estimate new value of xi
        x=ltitr(AKC,[K B-K*D],z);
        y0=x*C';
        if ~isempty(D),
            y0=y0+(D*z(:,ny+1:ny+nu).').';
        end
        psix0=ltitr(AKC.',C.',[1;zeros(Ncap,1)]);
        psix0=psix0(2:end,:);
        try
            xi(:,kexp) = pinv(psix0)*(z(:,1)-y0);
        catch
            warning('ident:estimation:X0EstFailed',...
                'Failed to estimate initial conditions. Check model stability. \nInitial state has been set to zero.')
            xi(:,kexp) = zeros(n,1);
        end
    end
end

%%*************************************************************************
function thpred = polypred(theta,m)
% Note that cross couplings between output
% channels are ignored in these calculations.

[ny,nu]= size(theta);
yna = pvget(theta,'OutputName');
una = pvget(theta,'InputName');
yu = pvget(theta,'OutputUnit');
uu = pvget(theta,'InputUnit');
for ky = 1:ny
    [a,b,c,d,f]=polydata(theta(ky,:));
    if nu>0
        ff=1;
        for ku=1:nu,
            bt=b(ku,:);
            for kku=1:nu,if kku~=ku,bt=conv(bt,f(kku,:));end,end
            bb(ku,:)=bt;
            ff=conv(ff,f(ku,:));
        end
        a=conv(conv(a,ff),d);c=conv(c,ff);
    else
        a=conv(a,d);
    end

    na=length(a);nc=length(c);nn=max(na,nc);
    a=[a,zeros(1,nn-na)];c=[c,zeros(1,nn-nc)];
    [f,g]=deconv(conv([1 zeros(1,m-1)],c),a);
    ng=length(g);
    if nu>0,
        df=conv(d,f);

        for ku=1:nu
            bf(ku,:)=conv(bb(ku,:),df);
        end
        nbf=length(bf(1,:));nn=max(ng,nbf);
        gg=[[g,zeros(1,nn-ng)];[bf,zeros(nu,nn-nbf)]];
    else
        gg=g;
    end
    th1 = idpoly(c,gg);
    th1 = pvset(th1,'InputName',[yna(ky);una],'OutputName',[yna{ky},'p'],...
        'InputUnit',[yu(ky);uu],'OutputUnit',yu(ky),'Ts',pvget(theta,'Ts'),...
        'TimeUnit',pvget(theta,'TimeUnit'),'EstimationInfo',pvget(theta,'EstimationInfo'));
    thpred{ky} =th1;
end
%%%************************************************************************
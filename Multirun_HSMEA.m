function Multirun_HSMEA
filename = mfilename;
motherfolderpath = which(filename);
motherfolder = fileparts(motherfolderpath);
cd(motherfolder);
Parameters;
addpath(genpath(motherfolder));
datafolder = [motherfolder,filesep,'Data'];
if ~exist(datafolder,'dir')
    mkdir(datafolder);
end

Problems = {
    'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7',...
    'WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9',...
    'DTLZ_1','DTLZ_2','DTLZ_3','DTLZ_4',...
    'WFG_1','WFG_2','WFG_3','WFG_4','WFG_5','WFG_6','WFG_7','WFG_8','WFG_9',...
    'C1_DTLZ1','C2_DTLZ2','C3_DTLZ4','CSI','WRM','GAA'
    };

numruns  = Global.numruns;

count = 1; Runpaths = []; ParamsAllruns = [];

for p = 1:numel(Problems)
    P = Problems{p};
    
    switch P 
        case 'CSI'
            Objs = 3;
            Global.FEmax = 700;
        case 'WRM'
            Objs = 5;
            Global.FEmax = 700;
        case 'GAA'
            Objs = 10;
            Global.FEmax = 700;
        otherwise
            Objs = [3 4 6 8 10];
            Global.FEmax = 300;
    end
   
    for m = Objs
        Global.M = m;
        Global.problem = P;
        prob = feval(P,m);
        Global.lower = prob.range(:,1)';
        Global.upper = prob.range(:,2)';
        Global.D = LoadDim(m,P);
        for run = 1:numruns
            Global.run = run;
            datapath = [datafolder,filesep,'HSMEA_',P,'_M',num2str(m),'_',num2str(run),'.mat'];
            Global.path = datapath;
            Global.filename = ['HSMEA_',P,'_M',num2str(m),'_',num2str(run),'.mat'];
            if ~exist(datapath,'file')
                Runpaths{count,1} = datapath;
                ParamsAllruns{count,1} = Global;
                count = count + 1;
            end
        end
    end
end

if ~isempty(Runpaths)
    cd(datafolder);
    parfor r = 1:numel(Runpaths)
        Global = ParamsAllruns{r};
        HSMEA(Global);
    end
end
end

function D = LoadDim(m,P)
if (m == 3 || m == 4) && ~isempty(regexp(P,'WFG', 'once'))
    D = 10;
elseif (m == 6 || m == 8) && ~isempty(regexp(P,'WFG', 'once'))
    D = 9;
elseif m == 10 && ~isempty(regexp(P,'WFG', 'once'))
    D = 11;
elseif ~isempty(regexp(P,'MaF', 'once'))
    if strcmpi(P,'MaF7')
        D = m + 19;
    elseif strcmpi(P,'MaF8') || strcmpi(P,'MaF9')
        D = 2;
    elseif strcmpi(P,'MaF13')
        D = 5;
    elseif strcmpi(P,'MaF14') || strcmpi(P,'MaF15')
        D = 3 * m + 5;
    else
        D = 10;
    end
elseif ~isempty(regexp(P,'CSI', 'once'))
    D = 7;
elseif ~isempty(regexp(P,'WRM', 'once'))
    D = 3;
elseif ~isempty(regexp(P,'GAA', 'once'))
    D = 27;
else
    D = max(10,m);
end
end
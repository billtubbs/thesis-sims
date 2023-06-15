function mkf_observer_sfunc(block)
%% Simulate a process observer in Simulink
%
%

setup(block);

%end mkf_observer_AFMM_sfunc



function setup(block)
%%   Set up the basic characteristics of the S-function block such as:
%      - Input ports
%      - Output ports
%      - Dialog parameters
%      - Options
%   Required         : Yes
%   C MEX counterpart: mdlInitializeSizes

% Register number of ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 2;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Get observer struct
obs = block.DialogPrm(1).Data;

% Input 1: u(k)
block.InputPort(1).Dimensions = obs.nu;
block.InputPort(1).DatatypeID = 0;  % double
block.InputPort(1).Complexity = 'Real';
block.InputPort(1).DirectFeedthrough = false;
block.InputPort(1).SamplingMode = 'Sample';

% Input 2: y(k)
block.InputPort(2).Dimensions = obs.ny;
block.InputPort(2).DatatypeID = 0;  % double
block.InputPort(2).Complexity = 'Real';
block.InputPort(2).DirectFeedthrough = false;
block.InputPort(2).SamplingMode = 'Sample';

% Output 1: x_est(k+1);
block.OutputPort(1).Dimensions = obs.n;
block.OutputPort(1).DatatypeID = 0; % double
block.OutputPort(1).Complexity = 'Real';
block.OutputPort(1).SamplingMode = 'Sample';

% Output 2: y_est(k+1)
block.OutputPort(2).Dimensions = obs.ny;
block.OutputPort(2).DatatypeID = 0; % double
block.OutputPort(2).Complexity = 'Real';
block.OutputPort(2).SamplingMode = 'Sample';

% Register parameters (number of parameters in the dialog box)
block.NumDialogPrms = 1;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [obs.Ts 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

% -----------------------------------------------------------------
% The MATLAB S-function uses an internal registry for all
% block methods. You should register all relevant methods
% (optional and required) as illustrated below. You may choose
% any suitable name for the methods and implement these methods
% as local functions within the same file. See comments
% provided for each function for more information.
% -----------------------------------------------------------------

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
%block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
%block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup



function DoPostPropSetup(block)
%% PostPropagationSetup:
%   Functionality    : Setup work areas and state variables. Can
%                      also register run-time methods here
%   Required         : No
%   C MEX counterpart: mdlSetWorkWidths

% Get observer struct
obs = block.DialogPrm(1).Data;

switch obs.type

    case {'KF', 'KFSS', 'LB'}  % obs with only double vars

        % Make data vectors containing all variables
        vec_double = get_obs_vars_vecs(obs);

        % Set number of Dwork blocks
        block.NumDworks = 1;

        % All dynamic variable data
        block.Dwork(1).Name            = 'vars_double';
        block.Dwork(1).Dimensions      = size(vec_double, 2);
        block.Dwork(1).DatatypeID      = 0;      % double
        block.Dwork(1).Complexity      = 'Real'; % real
        block.Dwork(1).UsedAsDiscState = true;

    case {'MKF', 'MKF_RODD', 'MKF_AFMM'}  % obs with double and int16

        % Make data vectors containing all variables
        [vec_double, vec_int16] = get_obs_vars_vecs(obs);

        % Set number of Dwork blocks
        block.NumDworks = 2;

        % All dynamic variable data
        block.Dwork(1).Name            = 'vars_double';
        block.Dwork(1).Dimensions      = size(vec_double, 2);
        block.Dwork(1).DatatypeID      = 0;      % double
        block.Dwork(1).Complexity      = 'Real'; % real
        block.Dwork(1).UsedAsDiscState = true;

        % All dynamic variable data
        block.Dwork(2).Name            = 'vars_int16';
        block.Dwork(2).Dimensions      = size(vec_int16, 2);
        block.Dwork(2).DatatypeID      = 4;      % int16
        block.Dwork(2).Complexity      = 'Real'; % real
        block.Dwork(2).UsedAsDiscState = true;

end



%end PostPropagationSetup



function InitializeConditions(block)
%%   Functionality    : Called at the start of simulation and if it is 
%                      present in an enabled subsystem configured to reset 
%                      states, it will be called when the enabled subsystem
%                      restarts execution to reset the states.
%   Required         : No
%   C MEX counterpart: mdlInitializeConditions  

% Get observer struct
obs = block.DialogPrm(1).Data;

switch obs.type
    
    case {'KF', 'KFSS', 'LB'}

        % Get data vectors containing all variables
        vec_double = get_obs_vars_vecs(obs);

        % Initialize Dwork memory vectors
        block.Dwork(1).Data = vec_double;

        % For debugging
        dlmwrite(sprintf('test-%s-double.csv', obs.label), vec_double, 'delimiter', ',');

    case {'MKF', 'MKF_RODD', 'MKF_AFMM'}

        % Get data vectors containing all variables
        [vec_double, vec_int16] = get_obs_vars_vecs(obs);

        % Initialize Dwork memory vectors
        block.Dwork(1).Data = vec_double;
        block.Dwork(2).Data = vec_int16;

        % For debugging
        dlmwrite(sprintf('test-%s-double.csv', obs.label), vec_double, 'delimiter', ',');
        dlmwrite(sprintf('test-%s-int16.csv', obs.label), vec_int16, 'delimiter', ',');

end


%end InitializeConditions



%function Start(block)
%%   Functionality    : Called once at start of model execution. If you
%                      have states that should be initialized once, this 
%                      is the place to do it.
%   Required         : No
%   C MEX counterpart: mdlStart%end Start

%end Start



function Outputs(block)
%%   Functionality    : Called to generate block outputs in
%                      simulation step
%   Required         : Yes
%   C MEX counterpart: mdlOutputs

% Get observer struct
obs = block.DialogPrm(1).Data;

switch obs.type

    case {'KF', 'KFSS', 'LB'}

        % Get variables data from Dwork memory
        vec_double = block.Dwork(1).Data;
        obs = set_obs_vars_vecs(obs, vec_double);

    case {'MKF', 'MKF_RODD', 'MKF_AFMM'}

        % Get variables data from Dwork memory
        vec_double = block.Dwork(1).Data;
        vec_int16 = block.Dwork(2).Data;
        obs = set_obs_vars_vecs(obs, vec_double, vec_int16);

    otherwise
        error('Value error: observer type not recognized')

end

% Output y_est(k+1)
block.OutputPort(1).Data = obs.xkp1_est;

% Output y_est(k+1)
block.OutputPort(2).Data = obs.ykp1_est;

%end Outputs



function Update(block)
%%   Functionality    : Called to update discrete states
%                      during simulation step
%   Required         : No
%   C MEX counterpart: mdlUpdate

% Get observer struct
obs = block.DialogPrm(1).Data;

% Inputs
uk = block.InputPort(1).Data;
yk = block.InputPort(2).Data;

% Check size of input vectors
assert(isequal(size(uk), [obs.nu 1]))
assert(isequal(size(yk), [obs.ny 1]))

switch obs.type

    case {'KF', 'KFSS', 'LB'}

        % Get variables data from Dwork memory
        vec_double = block.Dwork(1).Data;
        obs = set_obs_vars_vecs(obs, vec_double);

        % Update observer states
        obs = update_KF(obs, uk, yk);

        % Make data vectors containing all variables
        vec_double = get_obs_vars_vecs(obs);

        % For debugging
        dlmwrite(sprintf('test-%s-double.csv', obs.label), ...
            vec_double, 'delimiter', ',', '-append');

        % Update Dwork memory vectors
        block.Dwork(1).Data = vec_double;

    case {'MKF', 'MKF_RODD', 'MKF_AFMM'}

        % Get variables data from Dwork memory
        vec_double = block.Dwork(1).Data;
        vec_int16 = block.Dwork(2).Data;
        obs = set_obs_vars_vecs(obs, vec_double, vec_int16);

        % Update observer states
        switch obs.type
            case {'MKF', 'MKF_RODD'}
                  % For debugging
%                 if round(obs.filters{2}.P(2, 2), 4) == 1.3333
%                     disp('stop')
%                 end
                obs = update_MKF(obs, uk, yk);
            case 'MKF_AFMM'
                obs = update_AFMM(obs, uk, yk);
        end

        % Make data vectors containing all variables
        [vec_double, vec_int16] = get_obs_vars_vecs(obs);

        % For debugging
        dlmwrite(sprintf('test-%s-double.csv', obs.label), ...
            vec_double, 'delimiter', ',', '-append');
        dlmwrite(sprintf('test-%s-int16.csv', obs.label), ...
            vec_int16, 'delimiter', ',', '-append');

        % Update Dwork memory vectors
        block.Dwork(1).Data = vec_double;
        block.Dwork(2).Data = vec_int16;

    otherwise
        error('Value error: observer type not recognized')

end

%end Update



%function Derivatives(block)
%%   Functionality    : Called to update derivatives of
%                      continuous states during simulation step
%   Required         : No
%   C MEX counterpart: mdlDerivatives

%end Derivatives



function Terminate(block)
%%   Functionality    : Called at the end of simulation for cleanup
%   Required         : Yes
%   C MEX counterpart: mdlTerminate

%end Terminate



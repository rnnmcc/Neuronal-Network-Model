clear; clf; rng('default'); rng(17328461);

%% 

% Define network parameters: 

num_e=200; num_i=50;       % number of excitatory, inhibitory neurons
t_final=2000;               % Time (in ms) simulated.
sigma_e=0.05;              % standard deviation of excitatory input
% i_ext_e_hat = 1;           % average ecitatory input (mV)

sigma_i=0.05;               % standard deviation of inhibitory input
i_ext_i_hat = 0.7;            % average inhibitory input
i_ext_i=i_ext_i_hat*ones(num_i,1).*(1+sigma_i*randn(num_i,1));    

% i_ext_e_hat as a function of time
tstep1_i=200; tstep2_i= 350; step1_i=3.1; step2_i=3.1; step3_i = 3.1;

i_ext_e_hat = @(t) step1_i + (step2_i-step1_i)*(t>=tstep1_i) + (step3_i-step2_i)*(t>=tstep2_i);

i_ext_e = @(t) i_ext_e_hat(t)*ones(num_e,1).*(1+sigma_e*randn(num_e,1));     % store drive for each cell in vector


% M-Current- starts at step1_m before falling to step2_m, leading to
% hyperexcitability

step1_m = 1.2; step2_m = 0.1; tstep1_m = 1000;
g_m = @(t) step1_m + (step2_m-step1_m)*(t>=tstep1_m);

%% Synapses

g_hat_ee=0.1; g_hat_ei=0.5; g_hat_ie=0.5; g_hat_ii=0.5; % conductances for synapses (E to E etc.)
p_ee=0.5; p_ei=0.5; p_ie=0.5; p_ii=0.5;                   % connectivity of system
% p_ee=0.05; p_ei=0.05; p_ie=0.05; p_ii=0.05;

u_ee=rand(num_e,num_e); u_ei=rand(num_e,num_i);     % generating matrices of randoms in (0,1)
u_ie=rand(num_i,num_e); u_ii=rand(num_i,num_i);

g_ee=g_hat_ee*(u_ee<p_ee)/(num_e*p_ee); g_ei=g_hat_ei*(u_ei<p_ei)/(num_e*p_ei);     % pick out the elements less than p (should be approx p*num synapses) and multiply by g_hat/(num*p)
g_ie=g_hat_ie*(u_ie<p_ie)/(num_i*p_ie); g_ii=g_hat_ii*(u_ii<p_ii)/(num_i*p_ii);     
        % Consider, for example, the i-th e-cell and the j-th i-cell. the
        % probability that there is a synaptic connection at all from the
        % i-th e-cell to the j-th i-cell is p_ei. if there is such a
        % connection, its strength is g_hat_ei/(num_e*p_ei). Note that
        % num_e*p_ei is the expected number of excitatory inputs into an
        % inhibitory cell. Therefore dividing by this quantity has the
        % effect that the expected value of the total excitatory
        % conductance affecting an inhibitory cell is g_hat_ei. 
        
%% Gap Junctions

k_hat_ee=0.05; k_hat_ei=0.05; k_hat_ie=0.05; k_hat_ii=0.05;   % gap junction coupling strength (E to E etc.)
% p_gap_ee=0.025; p_gap_ei=0.0; p_gap_ie=0.0; p_gap_ii=0.0;    

p_gap_ee=0.025; p_gap_ei=0.05; p_gap_ie=0.05; p_gap_ii=0.025;

% GAP JUNCTION - II
u2_ee=rand(num_e,num_e); u2_ei=rand(num_e,num_i);     % generating matrices of randoms in (0,1)
u2_ie=rand(num_i,num_e); u2_ii=rand(num_i,num_i);

G_gap_ee=(u2_ee<p_gap_ee)/(num_e*p_gap_ee); G_gap_ei=(u2_ei<p_gap_ei)/(num_e*p_gap_ei);     % pick out the elements less than p (should be approx p*num synapses) and multiply by g_hat/(num*p)
G_gap_ie=(u2_ie<p_gap_ie)/(num_i*p_gap_ie); G_gap_ii=(u2_ii<p_gap_ii)/(num_i*p_gap_ii);

num_gap_ee = sum(G_gap_ee~=0,1)/(num_e*p_gap_ee);
num_gap_ei = sum(G_gap_ei~=0,1)/(num_e*p_gap_ei);
num_gap_ie = sum(G_gap_ie~=0,1)/(num_i*p_gap_ie);
num_gap_ii = sum(G_gap_ii~=0,1)/(num_i*p_gap_ii);



if p_gap_ii == 0
    G_gap_ii(isnan(G_gap_ii))=0;
    num_gap_ii(isinf(num_gap_ii))=0;
end
if p_gap_ee == 0
    G_gap_ee(isnan(G_gap_ee))=0;
    num_gap_ee(isinf(num_gap_ee))=0;
end
if p_gap_ei == 0
    G_gap_ei(isnan(G_gap_ei))=0;
    num_gap_ei(isinf(num_gap_ei))=0;
end
if p_gap_ie == 0
    G_gap_ie(isnan(G_gap_ie))=0;
    num_gap_ie(isinf(num_gap_ie))=0;
end


%% 

% tau_r ~ synaptic time rise constant, increasing leads to weaker synapse;
% tau_d ~ synaptic decay time constant; 
% tau_peak ~ time to peak
v_rev_e=0; v_rev_i=-75;         % reversal potentials
tau_r_e=0.5; tau_peak_e=0.5; %tau_d_e=3;         
tau_r_i=0.5; tau_peak_i=0.3; tau_d_i=9; 
 
dt=0.01;       % Time step used in solving the differential equations.
m_steps=round(t_final/dt);

times = linspace(0,t_final,m_steps);

% Introduce time dependence to taus - step at t = tstep, from step1 to
% step2

tstep = 200; step1 = 3.0; step2 = 3.0;
tau_d_e = @(t) step1 + (step2-step1)*(t>=tstep); % - (step2-step1)*(t>=200);

        
%% HH

% tau_dq_e=tau_d_q_function(tau_d_e(t),tau_r_e,tau_peak_e);
% tau_dq_i=tau_d_q_function(tau_d_i,tau_r_i,tau_peak_i);
dt05=dt/2; m_steps=round(t_final/dt);

% initialize dynamic variables
t0 = 0;

iv=rtm_init_w(i_ext_e(t0),rand(num_e,1),g_m(t0));
v_e=iv(:,1); m_e=m_e_inf(v_e); h_e=iv(:,2); n_e=iv(:,3); w=iv(:,4);
q_e=zeros(num_e,1); s_e=zeros(num_e,1);

v_i=-75*ones(num_i,1); m_i=m_i_inf(v_i); h_i=h_i_inf(v_i); n_i=n_i_inf(v_i); % Start all i-cells at -75 mv
z=zeros(num_i,1); q_i=z; s_i=z; 

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_e=0; t_e_spikes=[]; i_e_spikes=[];
num_spikes_i=0; t_i_spikes=[]; i_i_spikes=[];
lfp=zeros(m_steps+1,1); lfp(1)=mean(v_e);
lfp_i=zeros(m_steps+1,1); lfp(1)=mean(v_i);

v_e_store = zeros(num_e,m_steps);
v_i_store = zeros(num_i,m_steps);


%%
for k=1:m_steps
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
    tau_dq_e=tau_d_q_function(tau_d_e(t_old),tau_r_e,tau_peak_e);
    tau_dq_i=tau_d_q_function(tau_d_i,tau_r_i,tau_peak_i);

	v_e_inc=0.1*(-67-v_e)+80*n_e.^4.*(-100-v_e)+100*m_e.^3.*h_e.*(50-v_e) ...
               +(g_ee'*s_e).*(v_rev_e-v_e)+(g_ie'*s_i).*(v_rev_i-v_e)...
               +k_hat_ee*(G_gap_ee'*v_e-num_gap_ee'.*v_e)+k_hat_ie*(G_gap_ie'*v_i-num_gap_ie'.*v_e)...
               +g_m(t_old)*w.*(-100-v_e)+i_ext_e(t_old);

	n_e_inc=(n_e_inf(v_e)-n_e)./tau_n_e(v_e);
    h_e_inc=(h_e_inf(v_e)-h_e)./tau_h_e(v_e);
    q_e_inc=(1+tanh(v_e/10))/2.*(1-q_e)/0.1-q_e./tau_dq_e;
    s_e_inc=q_e.*(1-s_e)./tau_r_e-s_e./tau_d_e(t_old);
    w_inc=(w_inf(v_e)-w)./tau_w(v_e);
    
	v_i_inc=0.1*(-65-v_i)+9*n_i.^4.*(-90-v_i)+35*m_i.^3.*h_i.*(55-v_i) ...
               +(g_ei'*s_e).*(v_rev_e-v_i)+(g_ii'*s_i).*(v_rev_i-v_i) ...
               +k_hat_ei*(G_gap_ei'*v_e-num_gap_ei'.*v_i)+k_hat_ii*(G_gap_ii'*v_i-num_gap_ii'.*v_i) ... % NEW - gap junctions
               +i_ext_i;
           
	n_i_inc=(n_i_inf(v_i)-n_i)./tau_n_i(v_i);
    h_i_inc=(h_i_inf(v_i)-h_i)./tau_h_i(v_i);
    q_i_inc=(1+tanh(v_i/10))/2.*(1-q_i)/0.1-q_i./tau_dq_i;
    s_i_inc=q_i.*(1-s_i)./tau_r_i-s_i./tau_d_i;

	v_e_tmp=v_e+dt05*v_e_inc;
	n_e_tmp=n_e+dt05*n_e_inc;
	m_e_tmp=m_e_inf(v_e_tmp);
    h_e_tmp=h_e+dt05*h_e_inc;
    q_e_tmp=q_e+dt05*q_e_inc;   
    s_e_tmp=s_e+dt05*s_e_inc;
    w_tmp=w+dt05*w_inc;
	v_i_tmp=v_i+dt05*v_i_inc;
	n_i_tmp=n_i+dt05*n_i_inc;
	m_i_tmp=m_i_inf(v_i_tmp);
    h_i_tmp=h_i+dt05*h_i_inc;
    q_i_tmp=q_i+dt05*q_i_inc;   
    s_i_tmp=s_i+dt05*s_i_inc;    

	v_e_inc=0.1*(-67-v_e_tmp)+80*n_e_tmp.^4.*(-100-v_e_tmp)+100*m_e_tmp.^3.*h_e_tmp.*(50-v_e_tmp) ...
               +(g_ee'*s_e_tmp).*(v_rev_e-v_e_tmp)+(g_ie'*s_i_tmp).*(v_rev_i-v_e_tmp) ...
               +k_hat_ee*(G_gap_ee'*v_e_tmp-num_gap_ee'.*v_e_tmp)+k_hat_ie*(G_gap_ie'*v_i_tmp-num_gap_ie'.*v_e_tmp)...
               +g_m(t_old)*w_tmp.*(-100-v_e_tmp)+i_ext_e(t_old);     % +G_gap_ie*v_e_tmp-c_ie.*v_e_tmp;
	n_e_inc=(n_e_inf(v_e_tmp)-n_e_tmp)./tau_n_e(v_e_tmp);
    h_e_inc=(h_e_inf(v_e_tmp)-h_e_tmp)./tau_h_e(v_e_tmp);
    q_e_inc=(1+tanh(v_e_tmp/10))/2.*(1-q_e_tmp)/0.1-q_e_tmp./tau_dq_e;
    s_e_inc=q_e_tmp.*(1-s_e_tmp)./tau_r_e-s_e_tmp./tau_d_e(t_old);
    w_inc=(w_inf(v_e_tmp)-w_tmp)./tau_w(v_e_tmp);
    
	v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp)+35*m_i_tmp.^3.*h_i_tmp.*(55-v_i_tmp) ...
               +(g_ei'*s_e_tmp).*(v_rev_e-v_i_tmp)+(g_ii'*s_i_tmp).*(v_rev_i-v_i_tmp) ...
               +k_hat_ei*(G_gap_ei'*v_e_tmp-num_gap_ei'.*v_i_tmp)+k_hat_ii*(G_gap_ii'*v_i_tmp-num_gap_ii'.*v_i_tmp)...
               +i_ext_i;
	n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
    h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
    q_i_inc=(1+tanh(v_i_tmp/10))/2.*(1-q_i_tmp)/0.1-q_i_tmp./tau_dq_i;
    s_i_inc=q_i_tmp.*(1-s_i_tmp)./tau_r_i-s_i_tmp./tau_d_i;
        
    v_e_old=v_e;
    v_i_old=v_i;

	v_e=v_e+dt*v_e_inc;
	m_e=m_e_inf(v_e); h_e=h_e+dt*h_e_inc; n_e=n_e+dt*n_e_inc; 
    q_e=q_e+dt*q_e_inc;
    s_e=s_e+dt*s_e_inc;
	v_i=v_i+dt*v_i_inc;
    w=w+dt*w_inc;
    m_i=m_i_inf(v_i); h_i=h_i+dt*h_i_inc; n_i=n_i+dt*n_i_inc; 
    q_i=q_i+dt*q_i_inc;
    s_i=s_i+dt*s_i_inc;
    
    v_e_store(:,k) = v_e;   % stores for each individual neuron
    v_i_store(:,k) = v_i;
    
	% Determine which and how many e- and i-cells spiked in the current 
    % time step:

    which_e=find(v_e_old>-20 & v_e <=-20); which_i=find(v_i_old>-20 & v_i <=-20);
    l_e=length(which_e); l_i=length(which_i);
    if l_e>0 
        range=num_spikes_e+1:num_spikes_e+l_e; i_e_spikes(range)=which_e; 
        t_e_spikes(range)= ...
            ((-20-v_e(which_e))*(k-1)*dt+(v_e_old(which_e)+20)*k*dt)./ ...
                (-v_e(which_e)+v_e_old(which_e));
        num_spikes_e=num_spikes_e+l_e;
    end 
    if l_i>0 
        range=num_spikes_i+1:num_spikes_i+l_i; i_i_spikes(range)=which_i; 
        t_i_spikes(range)= ...
            ((-20-v_i(which_i))*(k-1)*dt+(v_i_old(which_i)+20)*k*dt)./ ...
                (-v_i(which_i)+v_i_old(which_i));
        num_spikes_i=num_spikes_i+l_i;
    end 
    
    lfp(k+1)=mean(v_e);
    lfp_i(k+1)=mean(v_i);

end


%% 

% plot the spike rastergram and mean voltage of E cells

% rastergram;


% plot(times,v_e_store(2,:)) %times, v_e_store(30,:))
% axis([0 400 -100 60])
% title('Single E-Cell Voltage Trace','Fontsize',20)
% xlabel('Time (ms)')
% ylabel('V (mV)')
% legend('E','I')


figure1 = figure();

% pre = tau_d_e(tstep-1); post = tau_d_e(tstep+1);
% i_pre = round(mean(i_ext_e(tstep_i-1)),2); i_post = round(mean(i_ext_e(tstep_i+1)),2);

if num_spikes_i>0, plot(t_i_spikes,i_i_spikes,'.b'); hold on; end
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes+num_i,'.r');  hold on; end
plot([0,t_final],[num_i+1/2,num_i+1/2],'--k','Linewidth',1);
hold off;
set(gca,'Fontsize',16); 
set(gca,'Ytick',[num_i,num_e+num_i]);
title({['Spike times']})
%     : N_e=' num2str(num_e) ' N_i=' num2str(num_i) ]})
xlabel('Time (ms)'); ylabel('Neuron Number');
%      ['i_e_x_t_ _e: ' num2str(step1_i) ' to ' num2str(step2_i)]})
%      ['g_i_i = ' num2str(g_hat_gap_ii) ' g_e_e = ' num2str(g_hat_gap_ee)] 
%      ['tau_r_e = ' num2str(tau_r_e) ' tau_p_e_a_k_ _e = ' num2str(tau_peak_e) ' tau_d_e = ' num2str(pre) ' => ' num2str(post) ]
%      ['tau_r_i = ' num2str(tau_r_i) ' tau_p_e_a_k_ _i = ' num2str(tau_peak_i) ' tau_d_i = ' num2str(tau_d_i)]})

% axis([-10,t_final,50.8,53.2]);
axis([-10,t_final,0,num_e+num_i+1]);      % set axis limits

% filepath = 'C:\Users\Ronan\Documents\College\Stage 4\Project\Non spike\Model\Simulations\Figures\Tau varied';
% filename = ['tau_r_i = ' num2str(tau_r_i) '.jpg'];
% saveas(figure1, fullfile(filepath,filename));

% plot average lfp


figure2 = figure();

subplot(211)
plot((0:m_steps)*dt,lfp,'-k','Linewidth',0.5);
set(gca,'Fontsize',16);
xlabel('t (ms)','Fontsize',15);
ylabel('V (mV)','Fontsize',15);
title('Mean LFP: E-cells')
axis([0,t_final,-100,50]);

subplot(212)
plot((0:m_steps)*dt,lfp_i,'-k','Linewidth',0.5);
set(gca,'Fontsize',16);
xlabel('t (ms)','Fontsize',15);
ylabel('V (mV)','Fontsize',15);
title('Mean LFP: I-cells')
axis([0,t_final,-100,50]);
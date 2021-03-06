function mpc = case33sj
%CASE33BW  Power flow data for 33 bus distribution system from Baran & Wu
%    Please see CASEFORMAT for details on the case file format.
%
%    Data from ...
%       M. E. Baran and F. F. Wu, "Network reconfiguration in distribution
%       systems for loss reduction and load balancing," in IEEE Transactions
%       on Power Delivery, vol. 4, no. 2, pp. 1401-1407, Apr 1989.
%       doi: 10.1109/61.25627
%       URL: http://doi.org/10.1109/61.25627

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [  %% (Pd and Qd are specified in kW & kVAr here, converted to MW & MVAr below)
	1	3	0	0	0	0	1	1	0	12.66	1	1	1;
	2	1	100	60	0	0	1	1	0	12.66	1	1.1	0.9;
	3	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	4	1	120	80	0	0	1	1	0	12.66	1	1.1	0.9;
	5	1	60	30	0	0	1	1	0	12.66	1	1.1	0.9;
	6	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	7	1	200	100	0	0	1	1	0	12.66	1	1.1	0.9;
	8	1	200	100	0	0	1	1	0	12.66	1	1.1	0.9;
	9	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	10	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	11	1	45	30	0	0	1	1	0	12.66	1	1.1	0.9;
	12	1	60	35	0	0	1	1	0	12.66	1	1.1	0.9;
	13	1	60	35	0	0	1	1	0	12.66	1	1.1	0.9;
	14	1	120	80	0	0	1	1	0	12.66	1	1.1	0.9;
	15	1	60	10	0	0	1	1	0	12.66	1	1.1	0.9;
	16	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	17	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	18	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	19	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	20	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	21	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	22	1	90	40	0	0	1	1	0	12.66	1	1.1	0.9;
	23	1	90	50	0	0	1	1	0	12.66	1	1.1	0.9;
	24	1	420	200	0	0	1	1	0	12.66	1	1.1	0.9;
	25	1	420	200	0	0	1	1	0	12.66	1	1.1	0.9;
	26	1	60	25	0	0	1	1	0	12.66	1	1.1	0.9;
	27	1	60	25	0	0	1	1	0	12.66	1	1.1	0.9;
	28	1	60	20	0	0	1	1	0	12.66	1	1.1	0.9;
	29	1	120	70	0	0	1	1	0	12.66	1	1.1	0.9;
	30	1	200	600	0	0	1	1	0	12.66	1	1.1	0.9;
	31	1	150	70	0	0	1	1	0	12.66	1	1.1	0.9;
	32	1	210	100	0	0	1	1	0	12.66	1	1.1	0.9;
	33	1	60	40	0	0	1	1	0	12.66	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	0	0	10	-10	1.00	100	1	10	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [  %% (r and x specified in ohms here, converted to p.u. below)
	1	2	0.0922	0.0470	0	0	0	0	0	0	1	-360	360;
	2	3	0.4930	0.2511	0	0	0	0	0	0	1	-360	360;
	3	4	0.3660	0.1864	0	0	0	0	0	0	1	-360	360;
	4	5	0.3811	0.1941	0	0	0	0	0	0	1	-360	360;
	5	6	0.8190	0.7070	0	0	0	0	0	0	1	-360	360;
	6	7	0.1872	0.6188	0	0	0	0	0	0	1	-360	360;
	7	8	0.7114	0.2351	0	0	0	0	0	0	1	-360	360;
	8	9	1.0300	0.7400	0	0	0	0	0	0	1	-360	360;
	9	10	1.0440	0.7400	0	0	0	0	0	0	1	-360	360;
	10	11	0.1966	0.0650	0	0	0	0	0	0	1	-360	360;
	11	12	0.3744	0.1238	0	0	0	0	0	0	1	-360	360;
	12	13	1.4680	1.1550	0	0	0	0	0	0	1	-360	360;
	13	14	0.5416	0.7129	0	0	0	0	0	0	1	-360	360;
	14	15	0.5910	0.5260	0	0	0	0	0	0	1	-360	360;
	15	16	0.7463	0.5450	0	0	0	0	0	0	1	-360	360;
	16	17	1.2890	1.7210	0	0	0	0	0	0	1	-360	360;
	17	18	0.7320	0.5740	0	0	0	0	0	0	1	-360	360;
	2	19	0.1640	0.1565	0	0	0	0	0	0	1	-360	360;
	19	20	1.5042	1.3554	0	0	0	0	0	0	1	-360	360;
	20	21	0.4095	0.4784	0	0	0	0	0	0	1	-360	360;
	21	22	0.7089	0.9373	0	0	0	0	0	0	1	-360	360;
	3	23	0.4512	0.3083	0	0	0	0	0	0	1	-360	360;
	23	24	0.8980	0.7091	0	0	0	0	0	0	1	-360	360;
	24	25	0.8960	0.7011	0	0	0	0	0	0	1	-360	360;
	6	26	0.2030	0.1034	0	0	0	0	0	0	1	-360	360;
	26	27	0.2842	0.1447	0	0	0	0	0	0	1	-360	360;
	27	28	1.0590	0.9337	0	0	0	0	0	0	1	-360	360;
	28	29	0.8042	0.7006	0	0	0	0	0	0	1	-360	360;
	29	30	0.5075	0.2585	0	0	0	0	0	0	1	-360	360;
	30	31	0.9744	0.9630	0	0	0	0	0	0	1	-360	360;
	31	32	0.3105	0.3619	0	0	0	0	0	0	1	-360	360;
	32	33	0.3410	0.5302	0	0	0	0	0	0	1	-360	360;
];
mpc.XY_data = [1	141	297
2	220	297
3	296	297
4	370	297
5	447	297
6	522	297
7	593	297
8	662	297
9	730	297
10	800	297
11	871	297
12	941	297
13	1011	297
14	1081	297
15	1149	297
16	1220	297
17	1290	297
18	1360	297
19	254	342
20	341	342
21	429	342
22	518	342
23	347	220
24	428	220
25	513	220
26	551	256
27	622	256
28	690	256
29	763	256
30	833	256
31	903	256
32	975	256
33	1044	256];
%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0	20	0;
];

%% convert branch impedances from Ohms to p.u.
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus; 
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch; 
Vbase = mpc.bus(1, 10) * 1e3;      %% in Volts
Sbase = mpc.baseMVA * 1e6;              %% in VA
mpc.branch(:, [3 4]) = mpc.branch(:, [3 4]) / (Vbase^2 / Sbase);
%% convert loads from kW to MW
mpc.bus(:, [3, 4]) = mpc.bus(:, [3, 4]) / 1e3;

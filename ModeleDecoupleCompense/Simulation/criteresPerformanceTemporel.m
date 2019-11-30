function [poles_des, wn, zeta] = criteresPerformanceTemporel(Mmax, ts, tp, tr10_90, tr0_100)
% [poles_des, wn, zeta] = criteresPerformanceTemporel(Mmax, ts, tp, tr10_90, tr0_100)
    zeta = cosd(atand(-pi/(log(Mmax/100))));
    if ts ~= 0 wn_ts = 4/(zeta*ts); else wn_ts = 0; end;
    if tp ~= 0 wn_tp = pi/(tp*sqrt(1-zeta^2)); else wn_tp = 0; end;
    if tr10_90 ~= 0 wn_tr10_90 = (1 + 1.1*zeta + 1.4*zeta^2)/tr10_90; else wn_tr10_90 = 0; end;
    if tr0_100 ~= 0 wn_tr0_100 = (pi-acos(zeta))/(tr0_100*sqrt(1-zeta^2)); else wn_tr0_100 = 0; end;
    wn = max([wn_ts wn_tp wn_tr10_90 wn_tr0_100]);
    poles_des(1) = -zeta*wn + wn*sqrt(1-zeta^2)*1i;
    poles_des(2) = -zeta*wn - wn*sqrt(1-zeta^2)*1i;
end


clear;
a = 2;
b = 9;
t = 0:0.1:30;
x = a * sin(0.5 * b * pi * t);
%x = a * sin(0.5 * b * t);
Ns = 2;
ts = t(1:Ns:end);
xs = x(1:Ns:end);
ln = [2 4 8];
e = zeros(length(ln));
var_th = zeros(1, length(ln));
var_eq = zeros(1, length(ln));
sqnr_th = zeros(1, length(ln));
sqnr_eq = zeros(1, length(ln));
decoded = zeros(length(ln), length(xs));
i = 1;
for l = ln
    [xq, v, dict, enc, eff, cr] = quantize(xs, l);
    err = abs(xq-xs);
    e(i) = mean(err);
    var_th(i) = v;
    var_eq(i) = var(err);
    sqnr_th(i) = 3 * (l^2);
    sqnr_eq(i) = max(x)^2 / var_eq(i);
    decoded(i,:) = huffmandeco(enc, dict);
    % disp(eff);
    % disp(cr);
    i = i + 1;
end
subplot(311);
plot(ln, e);
legend('mean error'); 
subplot(312)
plot(ln, var_th, ln, var_eq);
legend('variance th', 'variance eq');
subplot(313);
plot(ln, sqnr_th, ln, sqnr_eq);
legend('sqnr th','sqnr eq');

for i = 1:length(ln)
    subplot(length(ln), 1, i);
    plot(t, x, ts, decoded(i,:));
%    xlim([0, 10]);
    legend(sprintf('%d levels', 2^i));
end






function [xq, variance, dict, enc, eff, cr] =  quantize(xs, l)
    vmin = min(xs);
    vmax = max(xs);
    delta = (vmax - vmin) / l;
    levels = (vmin +(delta / 2)):delta:(vmax - (delta / 2));
    prob = zeros(1, l);
    xq = xs;
    i = 1;
    format long;
    for le = levels
        ind = find((i == 1 & xq == le - (delta/2)) | xq > le - (delta/2) & xq <= le + (delta/2));
        prob(i) = length(ind) / length(xs);
        xq(ind) = le;
        i = i + 1;
    end

    variance = delta * delta / 12;
    [dict,avgl] = huffmandict(levels, prob);
    enc = huffmanenco(xq, dict);
    lmin =  - sum(prob.* log2(prob));
    eff = lmin / avgl;
    cr = ceil(log2(l)) / avgl;
end

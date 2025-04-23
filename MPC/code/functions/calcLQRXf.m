function [Xf_H, Xf_h] = calcLQRXf(sysd, xlb, xub, ulb, uub, Q, R)
model = LTISystem(sysd);
model.x.min = xlb;
model.x.max = xub;

model.u.min = ulb;
model.u.max = uub;

model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);

Tset = model.LQRSet;

Xf_H = Tset.A;
Xf_h = Tset.b;

end
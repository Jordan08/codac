import py.codac4matlab.*

% [B-q3-beg]
M = {IntervalVector([1.5,2.5]),IntervalVector([3,1]), IntervalVector([2,2]), IntervalVector([2.5,3]), IntervalVector([3.5,2]), IntervalVector([4,1]), IntervalVector([1.5,0.5])};

for i = 1:numel(M)
    M{i}.inflate(0.05);
end

a1 = IntervalVector({{1.25,3.},{1.6,2.75}});
a2 = IntervalVector({{2,3.5},{0.6,1.2}});
a3 = IntervalVector({{1.1,3.25},{0.2,1.4}});

ctc_constell = MyCtc(M);

a1 = ctc_constell.contract(a1)
a2 = ctc_constell.contract(a2)
a3 = ctc_constell.contract(a3)
% [B-q3-end]

% [B-q4-beg]
x_truth = Vector([2,1,PI/6]);
% [B-q4-end]

% [B-q5-beg]
DefaultFigure().draw_tank(x_truth, 0.4, StyleProperties({Color().black(),Color().yellow()}))
for i = 1:numel(M)
  DefaultFigure().draw_box(M{i}, StyleProperties({Color().dark_green(),Color().green()}))
end

DefaultFigure().set_axes(axis(1,Interval([1.,4.5])), axis(2,Interval([0.,3.5])));
% [B-q5-end]
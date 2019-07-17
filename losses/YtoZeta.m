function [zeta] = YtoZeta()

  zeta = (1 - (1 + Y * (1 - p2/pt2))^((1 - gamma)/gamma)) / ((gamma - 1) / 2 * M2 ^ 2)

end

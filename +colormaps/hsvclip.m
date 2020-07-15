function cm = cm_hsvclip

cm = hsv(10);
cm = kron(cm,ones(15,1));
cm = [zeros(1,3); cm; ones(1,3)];
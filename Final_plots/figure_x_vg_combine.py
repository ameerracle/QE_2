# Recovered from IPython history - session cell ending line 140517
# Combine the |V|/G panels into one image using the laplacian-merge style
import matplotlib.image as mpimg
# Panel 1: 2x3 grid image, Panel 2: 1x3 grid image
vg_a = mpimg.imread(base_dir / 'Final_plots' / 'Figure_abs_V_over_G.png')
vg_b = mpimg.imread(base_dir / 'Final_plots' / 'Figure_abs_V_over_G_by_ads.png')

# Match styling used for the Laplacian combine (consistent column width, small v-space)
col_width = 18
col_height = 6
fig = plt.figure(figsize=(col_width, col_height * 2 + 4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.67, 1.0], hspace=0.04)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])

ax_a.imshow(vg_a, aspect='auto')
ax_a.axis('off')
ax_a.set_title('a)', loc='left', fontsize=18, fontweight='bold', pad=10)

ax_b.imshow(vg_b, aspect='auto')
ax_b.axis('off')
ax_b.set_title('b)', loc='left', fontsize=18, fontweight='bold', pad=10)

fig.tight_layout(pad=0.4)
fig.savefig(base_dir / 'Final_plots' / 'Figure_5.png', dpi=600, bbox_inches='tight')
#plt.close(fig)

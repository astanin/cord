#!/usr/bin/env ruby

# plot phi colormap and interface position

require 'Tioga/FigureMaker'

class CordPlots

	include Math
	include Tioga
	include FigureConstants
    
	def t
		@figure_maker
	end

	def initialize
		@figure_maker = FigureMaker.default
		t.def_eval_function { |str| eval(str) }

		@phi_range = [0.7, 0.85]
		@phi_colormap = t.create_gradient_colormap(
			'starting_H' => 330,
			'starting_S' => 0.24,
			'starting_L' => 0.9,
			'ending_H' => 358,
			'ending_S' => 0.45,
			'ending_L' => 0.5
		)

		@c_range = [0.7, 1.0]
		@c_colormap = t.create_gradient_colormap(
			'starting_H' => 15,
			'starting_S' => 1.0,
			'starting_L' => 0.4,
			'ending_H' => 270,
			'ending_S' => 0.3,
			'ending_L' => 0.8
		)

		@loaded_data = false

		@data_filename = ARGV[0]
		t.auto_refresh_filename = @data_filename
		if ARGV.length == 2
			t.save_dir = ARGV[1]
		else
			t.save_dir = 'figures_out'
		end

		t.auto_refresh_filename=true
		t.add_model_number = false
		#t.tex_preview_paper_width="18cm"
		#t.tex_preview_paper_height="15cm"
		#t.tex_preview_figure_width="18cm"
		#t.tex_preview_figure_height="15cm"
		t.tex_preview_preamble += "\n\\include{color_names}\n"
		get_data unless @loaded_data
		c_filename="OxygenColormap_%07d" % (@time*10).floor
		phi_filename="PhiColormap_%07d" % (@time*10).floor
		@c_plot=t.def_figure(c_filename) {
			plot_c_colormap
		}
		@phi_plot=t.def_figure(phi_filename) {
			plot_phi_colormap
		}
		t.def_enter_page_function { enter_page }
	end

	def enter_page
		t.page_setup(5*72,5*72)
		t.set_frame_sides(0.15,0.9,0.90,0.10)
		t.set_default_font_size(14)
	end

	def choose_margins
		@aspect=(@xvar.max-@xvar.min)/(@yvar.max-@yvar.min);
		if (@aspect >= 1.0)
			@ml=0.0
			@mr=0.3
			@mlr=0.0
			@mtb=0.5*(1-(1-@ml-@mr)/@aspect)
		else
			@ml=0.0
			@mr=0.3
			@mtb=0.15
			@mlr=0.5*(1-@ml-@mr-@aspect*(1-2*@mtb))
		end
	end

	def plot_phi_colormap
		get_data unless @loaded_data
		title="$t=%5.1f$" % @time
		choose_margins
		t.subplot('left_margin'=>@ml+@mlr, 'right_margin'=>@mr+@mlr,
			'top_margin'=>@mtb, 'bottom_margin'=>@mtb) {
			t.do_box_labels(title, '$x$', '$y$')
			show_colormap(@phi_image, @phi_colormap)
			show_isolines(@psi, [0])
		}
		t.subplot('left_margin'=>(1-0.3*(@ml+@mr)), 'top_margin'=>@mtb,
			'bottom_margin'=>@mtb) {
			show_colorbar("$\\phi$", @phi_range, @phi_colormap)
		}
	end

	def plot_c_colormap
		get_data unless @loaded_data
		title="$t=%5.1f$" % @time
		choose_margins
		t.subplot('left_margin'=>@ml+@mlr, 'right_margin'=>@mr+@mlr,
			'top_margin'=>@mtb, 'bottom_margin'=>@mtb) {
			t.do_box_labels(title, '$x$', '$y$')
			show_colormap(@c_image, @c_colormap)
			show_isolines(@psi, [0])
		}
		t.subplot('left_margin'=>(1-0.3*(@ml+@mr)), 'top_margin'=>@mtb,
			'bottom_margin'=>@mtb) {
			show_colorbar("$c$", @c_range, @c_colormap)
		}
	end

	def show_colormap(var_image, colormap)
		t.line_width *= 0.5
		t.xaxis_line_width *= 0.5
		t.yaxis_line_width *= 0.5
		t.xaxis_major_tick_width *= 0.5
		t.xaxis_minor_tick_width *= 0.5
		t.yaxis_major_tick_width *= 0.5
		t.yaxis_minor_tick_width *= 0.5
		t.show_plot([@xvar.min,@xvar.max, @yvar.max,@yvar.min]) do
			t.fill_color=Wheat
			t.fill_frame
			t.show_image('data'=>var_image,
				'w' => @xdim, 'h' => @ydim,
				'll' => [ @xvar.min, @yvar.max ],
				'lr' => [ @xvar.max, @yvar.max ],
				'ul' => [ @xvar.min, @yvar.min ],
				'colormap' => colormap,
				'interpolate' => false)
		end
	end

	def show_isolines(var, levels)
		t.stroke_color=Black
		t.line_width=1
		dest_xs = Dvector.new; dest_ys = Dvector.new; gaps = Array.new
		dict={ 'dest_xs' => dest_xs, 'dest_ys' => dest_ys,
			'gaps' => gaps,
			'data' => var, 'xs' => @xvar, 'ys' => @yvar
#			/* , 'method' => 'conrec' */
			}
		levels.each { |level|
			dict['level'] = level
			t.make_contour(dict)
			t.append_points_with_gaps_to_path(
				dest_xs, dest_ys, gaps, false)
			t.stroke
		}
	end

	def show_colorbar(title, range, colormap)
		t.line_width *= 0.5
		t.xaxis_line_width *= 0.5
		t.yaxis_line_width *= 0.5
		t.xaxis_major_tick_width *= 0.5
		t.xaxis_minor_tick_width *= 0.5
		t.yaxis_major_tick_width *= 0.5
		t.yaxis_minor_tick_width *= 0.5
		xmin=0.0; xmid=0.5; xmax=1.0;
		t.xaxis_type = AXIS_LINE_ONLY
		t.xaxis_loc = BOTTOM
		t.top_edge_type = AXIS_LINE_ONLY
		t.yaxis_loc = RIGHT
		t.ylabel_side = BOTTOM
		t.yaxis_type = AXIS_WITH_TICKS_AND_NUMERIC_LABELS
		t.left_edge_type = AXIS_WITH_TICKS_ONLY
		t.ylabel_shift -= 0.1
		t.yaxis_major_tick_length *= 0.6
		t.yaxis_minor_tick_length *= 0.5
		t.show_title(title)
		#t.show_ylabel(title)
		t.no_ylabel
		t.show_plot([xmin,xmax,range[1],range[0]]) do
			t.axial_shading('start_point' => [xmid,range[0]],
				'end_point' => [xmid,range[1]],
				'colormap' => colormap)
		end
	end

	def get_data
		read_data
		t.model_number=(@time*10.0).floor
		@phi_image=t.create_image_data(@phi,
			'min_value' => @phi_range[0],
			'max_value' => @phi_range[1])
		@c_image=t.create_image_data(@c,
			'min_value' => @c_range[0],
			'max_value' => @c_range[1])
		@loaded_data = true
	end

	def read_data
		f=File.new(@data_filename,"r")
		l=f.readline
		@time=l.split()[1].to_f
		puts "DEBUG: time=#{@time}\n"
		l=f.readline
		@xdim=l.split()[1].to_i
		@ydim=l.split()[2].to_i
		puts "DEBUG: dimensions=#{@xdim}x#{@ydim}\n"
		@xvar=Dvector.new(0)
		@yvar=Dvector.new(0)
		x=f.readline.split.slice(1..-1)
		i=0; x.each{ |e| @xvar[i]=e.to_f; i=i.next }
		y=f.readline.split.slice(1..-1)
		i=0; y.each{ |e| @yvar[i]=e.to_f; i=i.next }
		puts "DEBUG: bounds=[#{@xvar.min}, #{@xvar.max}, " +
			"#{@yvar.min}, #{@yvar.max}]\n"
		@c=Dtable.new(@xdim,@ydim)
		@phi=Dtable.new(@xdim,@ydim)
		@psi=Dtable.new(@xdim,@ydim)
		# read variable names and positions
		vars=f.readline.split.collect { |v| v.sub('#','') }
		idxs=Hash.new
		vars.each { |v| idxs[v]=vars.index(v) }
		# read data
		k=0
		f.each_line do |l| 
			l=l.strip
			# skip comments and empty lines
			if (l[0..0] != '#') && (l !~ /^[ \t]*$/)
				i=k.divmod(@ydim)[0]
				j=k.divmod(@ydim)[1]
				a=l.split
				@c[j,i]  =a[idxs["c"]].to_f
				@phi[j,i]=a[idxs["phi"]].to_f
				@psi[j,i]=a[idxs["psi"]].to_f
				k=k.next
			end
		end
		f.close
	end
end

class CordPlotsBatch < CordPlots

	def run
		0.upto(t.num_figures-1) do |i|
			t.make_preview_pdf(i)
		end
	end

end

CordPlotsBatch.new.run


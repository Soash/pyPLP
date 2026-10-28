from pyPLP.plp_generation.backbone_optimizer import optimize_backbone_candidate
from pyPLP.plp_generation.random_backbone import generate_random_backbone
from pyPLP.plp_generation.segmented_backbone import build_initial_segmented_backbone

shared_counter = None
counter_lock = None
shared_stop_event = None


def init_worker(counter, lock, stop_event=None):
    """Initializes global shared counter, lock, and optional stop event for worker process."""
    global shared_counter, counter_lock, shared_stop_event
    shared_counter = counter
    counter_lock = lock
    shared_stop_event = stop_event


def worker_screen_batch(args):
    """Worker task that stochastically designs and tests random backbones until
    an open PLP is found or stop_event is set.
    """
    global shared_counter, counter_lock, shared_stop_event
    (
        left_arm,
        right_arm,
        backbone_len,
        max_homo,
        max_hairpin_tm,
        min_plp_mfe_dg,
        max_paired_arm_bases,
        max_paired_total_bases,
    ) = args[:8]
    arch_cfg = args[8] if len(args) > 8 else None

    local_ticks = 0

    while True:
        if shared_stop_event is not None and shared_stop_event.is_set():
            if local_ticks > 0 and counter_lock is not None and shared_counter is not None:
                with counter_lock:
                    shared_counter.value += local_ticks
            return None

        local_ticks += 1
        if local_ticks >= 25:
            if counter_lock is not None and shared_counter is not None:
                with counter_lock:
                    shared_counter.value += local_ticks
            local_ticks = 0

        if arch_cfg:
            initial_bb, prot_indices = build_initial_segmented_backbone(
                spacer1_len=arch_cfg.get("spacer1_len", 4),
                rca_primer_binding_len=arch_cfg.get("rca_primer_binding_len", 20),
                spacer2_len=arch_cfg.get("spacer2_len", 4),
                probe_binding_len=arch_cfg.get("probe_binding_len", 20),
                spacer3_len=arch_cfg.get("spacer3_len", 0),
                rest_len=arch_cfg.get("rest_len", 4),
                restriction_site=arch_cfg.get("restriction_site", ""),
                universal_rca_primer_seq=arch_cfg.get("universal_rca_primer_seq", ""),
                universal_probe_seq=arch_cfg.get("universal_probe_seq", ""),
                min_tm=arch_cfg.get("min_rca_primer_tm", 40.0),
                max_homo=max_homo,
                max_hairpin_tm=max_hairpin_tm,
                re_mode=arch_cfg.get("re_mode", "standalone"),
                re_probe_offset=arch_cfg.get("re_probe_offset", 2),
                crispr_mode=bool(arch_cfg.get("crispr_mode", False)),
                rest_tail=arch_cfg.get("rest_tail", ""),
            )
        else:
            initial_bb = generate_random_backbone(backbone_len)
            prot_indices = None

        candidate = optimize_backbone_candidate(
            initial_bb,
            left_arm,
            right_arm,
            max_homo=max_homo,
            max_hairpin_tm=max_hairpin_tm,
            min_plp_mfe_dg=min_plp_mfe_dg,
            max_paired_arm_bases=max_paired_arm_bases,
            max_paired_total_bases=max_paired_total_bases,
            arch_cfg=arch_cfg,
            protected_bb_indices=prot_indices,
        )

        if candidate is not None:
            if local_ticks > 0 and counter_lock is not None and shared_counter is not None:
                with counter_lock:
                    shared_counter.value += local_ticks
            return candidate

-- ============================================================================
-- TITLE: LaLiga outfield player-match extraction and player-season construction
-- AUTHOR: Reconstructed from staged Wyscout schema and analytical pipeline
-- PURPOSE:
--   1) Restrict the staged Wyscout warehouse to historical LaLiga data
--   2) Join competition, season, player, and match-player statistics tables
--   3) Exclude goalkeeper rows using goalkeeper-specific activity markers
--   4) Aggregate match-player rows to a player-season analytical dataset
--   5) Construct interpretable and comparable football features for analysis
--
-- OUTPUT GRAIN:
--   One row per player_id x season_id x competition_id
--
-- NOTES:
--   - This query follows the logic used in the R pipeline that later computes
--     distance-based similarity structures.
--   - The LaLiga restriction is treated as a design choice, not a limitation.
--   - The goalkeeper exclusion is performed behaviorally through goalkeeper
--     metrics, not solely through role labels.
--   - A minimum exposure threshold of 600 minutes is enforced to reduce noise.
-- ============================================================================

with

/* ---------------------------------------------------------------------------
   1) Identify the seasons belonging to the target competition universe
   --------------------------------------------------------------------------- */
target_competitions as (
    select
        c.competition_id,
        c.competition_name,
        c.competition_area_id,
        c.competition_area_alpha_2_code,
        c.competition_area_alpha_3_code,
        c.competition_area_name,
        c.competition_format,
        c.competition_type,
        c.competition_category,
        c.competition_gender,
        c.competition_division_level
    from prd.staging_wyscout.stg_wyscout__competitions as c
    where lower(c.competition_name) like '%laliga%'
      and lower(coalesce(c.competition_gender, 'male')) = 'male'
      and coalesce(c.competition_division_level, 1) = 1
),

target_seasons as (
    select
        s.season_id,
        s.season_name,
        s.season_start_date,
        s.season_end_date,
        s.season_is_active,
        s.competition_id
    from prd.staging_wyscout.stg_wyscout__seasons as s
    inner join target_competitions as tc
        on s.competition_id = tc.competition_id
),

/* ---------------------------------------------------------------------------
   2) Build the staged match-player base table at one row per player x match
   --------------------------------------------------------------------------- */
match_player_base as (
    select
        mps.player_id,
        p.player_name,
        p.player_short_name,
        p.player_name_unidecoded,
        p.player_height,
        p.player_weight,
        p.player_birth_date,
        p.player_birth_area_name,
        p.player_passport_area_name,
        p.role_name,
        p.role_code_2,
        p.role_code_3,
        p.player_prefered_foot,
        p.player_gender,
        p.player_status,

        mps.match_id,
        mps.competition_id,
        tc.competition_name,
        tc.competition_area_name,
        tc.competition_format,
        tc.competition_type,
        tc.competition_category,
        tc.competition_gender,
        tc.competition_division_level,

        mps.season_id,
        ts.season_name,
        ts.season_start_date,
        ts.season_end_date,
        ts.season_is_active,

        mps.round_id,

        -- exposure and usage
        mps.total_matches,
        mps.total_matches_in_start,
        mps.total_matches_substituted,
        mps.total_matches_coming_off,
        mps.total_minutes_on_field,
        mps.total_minutes_tagged,

        -- attacking output
        mps.total_goals,
        mps.total_assists,
        mps.total_shots,
        mps.total_head_shots,
        mps.total_shots_on_target,
        mps.total_shots_blocked,
        mps.total_penalties,
        mps.total_successful_penalties,
        mps.total_touch_in_box,
        mps.total_offsides,
        mps.total_xg_shot,
        mps.total_xg_assist,

        -- passing and creation
        mps.total_passes,
        mps.total_successful_passes,
        mps.total_smart_passes,
        mps.total_successful_smart_passes,
        mps.total_passes_to_final_third,
        mps.total_successful_passes_to_final_third,
        mps.total_crosses,
        mps.total_successful_crosses,
        mps.total_forward_passes,
        mps.total_successful_forward_passes,
        mps.total_back_passes,
        mps.total_successful_back_passes,
        mps.total_through_passes,
        mps.total_successful_through_passes,
        mps.total_key_passes,
        mps.total_successful_key_passes,
        mps.total_vertical_passes,
        mps.total_successful_vertical_passes,
        mps.total_long_passes,
        mps.total_successful_long_passes,
        mps.total_progressive_passes,
        mps.total_successful_progressive_passes,
        mps.total_received_pass,
        mps.total_shot_assists,
        mps.total_shot_on_target_assists,
        mps.total_second_assists,
        mps.total_third_assists,
        mps.total_linkup_plays,
        mps.total_successful_linkup_plays,

        -- carrying and ball progression
        mps.total_dribbles,
        mps.total_successful_dribbles,
        mps.total_progressive_run,
        mps.total_accelerations,

        -- defensive activity and pressing
        mps.total_defensive_actions,
        mps.total_successful_defensive_actions,
        mps.total_interceptions,
        mps.total_recoveries,
        mps.total_opponent_half_recoveries,
        mps.total_dangerous_opponent_half_recoveries,
        mps.total_counterpressing_recoveries,
        mps.total_pressing_duels,
        mps.total_pressing_duels_won,
        mps.total_clearances,
        mps.total_sliding_tackles,
        mps.total_successful_sliding_tackles,

        -- duels and physical contests
        mps.total_duels,
        mps.total_duels_won,
        mps.total_defensive_duels,
        mps.total_defensive_duels_won,
        mps.total_offensive_duels,
        mps.total_offensive_duels_won,
        mps.total_aerial_duels,
        mps.total_aerial_duels_won,
        mps.total_field_aerial_duels,
        mps.total_field_aerial_duels_won,
        mps.total_loose_ball_duels,
        mps.total_loose_ball_duels_won,

        -- discipline and ball security
        mps.total_fouls,
        mps.total_fouls_suffered,
        mps.total_yellow_cards,
        mps.total_red_cards,
        mps.total_direct_red_cards,
        mps.total_losses,
        mps.total_own_half_losses,
        mps.total_dangerous_own_half_losses,
        mps.total_missed_balls,

        -- goalkeeper activity markers retained only for exclusion logic
        mps.total_goal_kicks,
        mps.total_goal_kicks_short,
        mps.total_goal_kicks_long,
        mps.total_successful_goal_kicks,
        mps.total_gk_clean_sheets,
        mps.total_gk_conceded_goals,
        mps.total_gk_shots_against,
        mps.total_gk_exits,
        mps.total_gk_successful_exits,
        mps.total_gk_aerial_duels,
        mps.total_gk_aerial_duels_won,
        mps.total_gk_saves,
        mps.total_xg_save,

        mps._last_modified_at,
        mps._last_imported_at
    from prd.staging_wyscout.stg_wyscout__matches_players_stats as mps
    inner join target_seasons as ts
        on mps.season_id = ts.season_id
    inner join target_competitions as tc
        on mps.competition_id = tc.competition_id
    left join prd.staging_wyscout.stg_wyscout__players as p
        on mps.player_id = p.player_id
),

/* ---------------------------------------------------------------------------
   3) Exclude goalkeeper observations using goalkeeper-specific actions
   --------------------------------------------------------------------------- */
outfield_match_player_base as (
    select
        mpb.*
    from match_player_base as mpb
    where coalesce(mpb.total_gk_shots_against, 0) = 0
      and coalesce(mpb.total_gk_saves, 0) = 0
      and coalesce(mpb.total_gk_conceded_goals, 0) = 0
      and coalesce(mpb.total_gk_clean_sheets, 0) = 0
      and coalesce(mpb.total_gk_exits, 0) = 0
),

/* ---------------------------------------------------------------------------
   4) Aggregate match-player rows to player-season rows
   --------------------------------------------------------------------------- */
player_season_base as (
    select
        omp.player_id,
        max(omp.player_name) as player_name,
        max(omp.player_short_name) as player_short_name,
        max(omp.player_name_unidecoded) as player_name_unidecoded,
        max(omp.player_height) as player_height,
        max(omp.player_weight) as player_weight,
        max(omp.player_birth_date) as player_birth_date,
        max(omp.player_birth_area_name) as player_birth_area_name,
        max(omp.player_passport_area_name) as player_passport_area_name,
        max(omp.role_name) as role_name,
        max(omp.role_code_2) as role_code_2,
        max(omp.role_code_3) as role_code_3,
        max(omp.player_prefered_foot) as player_prefered_foot,

        omp.season_id,
        max(omp.season_name) as season_name,
        max(omp.season_start_date) as season_start_date,
        max(omp.season_end_date) as season_end_date,

        omp.competition_id,
        max(omp.competition_name) as competition_name,
        max(omp.competition_area_name) as competition_area_name,
        max(omp.competition_format) as competition_format,
        max(omp.competition_type) as competition_type,
        max(omp.competition_category) as competition_category,

        count(distinct omp.match_id) as matches,
        sum(coalesce(omp.total_matches_in_start, 0)) as starts,
        sum(coalesce(omp.total_matches_substituted, 0)) as matches_substituted,
        sum(coalesce(omp.total_matches_coming_off, 0)) as matches_coming_off,
        sum(coalesce(omp.total_minutes_on_field, 0)) as minutes,
        sum(coalesce(omp.total_minutes_tagged, 0)) as minutes_tagged,

        -- core totals
        sum(coalesce(omp.total_goals, 0)) as goals,
        sum(coalesce(omp.total_assists, 0)) as assists,
        sum(coalesce(omp.total_shots, 0)) as shots,
        sum(coalesce(omp.total_head_shots, 0)) as head_shots,
        sum(coalesce(omp.total_shots_on_target, 0)) as shots_on_target,
        sum(coalesce(omp.total_shots_blocked, 0)) as shots_blocked,
        sum(coalesce(omp.total_penalties, 0)) as penalties,
        sum(coalesce(omp.total_successful_penalties, 0)) as successful_penalties,
        sum(coalesce(omp.total_touch_in_box, 0)) as touches_box,
        sum(coalesce(omp.total_offsides, 0)) as offsides,
        sum(coalesce(omp.total_xg_shot, 0.0)) as xg_shot,
        sum(coalesce(omp.total_xg_assist, 0.0)) as xg_assist,

        -- passing and creation
        sum(coalesce(omp.total_passes, 0)) as passes,
        sum(coalesce(omp.total_successful_passes, 0)) as successful_passes,
        sum(coalesce(omp.total_smart_passes, 0)) as smart_passes,
        sum(coalesce(omp.total_successful_smart_passes, 0)) as successful_smart_passes,
        sum(coalesce(omp.total_passes_to_final_third, 0)) as passes_to_final_third,
        sum(coalesce(omp.total_successful_passes_to_final_third, 0)) as successful_passes_to_final_third,
        sum(coalesce(omp.total_crosses, 0)) as crosses,
        sum(coalesce(omp.total_successful_crosses, 0)) as successful_crosses,
        sum(coalesce(omp.total_forward_passes, 0)) as forward_passes,
        sum(coalesce(omp.total_successful_forward_passes, 0)) as successful_forward_passes,
        sum(coalesce(omp.total_back_passes, 0)) as back_passes,
        sum(coalesce(omp.total_successful_back_passes, 0)) as successful_back_passes,
        sum(coalesce(omp.total_through_passes, 0)) as through_passes,
        sum(coalesce(omp.total_successful_through_passes, 0)) as successful_through_passes,
        sum(coalesce(omp.total_key_passes, 0)) as key_passes,
        sum(coalesce(omp.total_successful_key_passes, 0)) as successful_key_passes,
        sum(coalesce(omp.total_vertical_passes, 0)) as vertical_passes,
        sum(coalesce(omp.total_successful_vertical_passes, 0)) as successful_vertical_passes,
        sum(coalesce(omp.total_long_passes, 0)) as long_passes,
        sum(coalesce(omp.total_successful_long_passes, 0)) as successful_long_passes,
        sum(coalesce(omp.total_progressive_passes, 0)) as progressive_passes,
        sum(coalesce(omp.total_successful_progressive_passes, 0)) as successful_progressive_passes,
        sum(coalesce(omp.total_received_pass, 0)) as received_passes,
        sum(coalesce(omp.total_shot_assists, 0)) as shot_assists,
        sum(coalesce(omp.total_shot_on_target_assists, 0)) as shot_on_target_assists,
        sum(coalesce(omp.total_second_assists, 0)) as second_assists,
        sum(coalesce(omp.total_third_assists, 0)) as third_assists,
        sum(coalesce(omp.total_linkup_plays, 0)) as linkup_plays,
        sum(coalesce(omp.total_successful_linkup_plays, 0)) as successful_linkup_plays,

        -- carrying and ball progression
        sum(coalesce(omp.total_dribbles, 0)) as dribbles,
        sum(coalesce(omp.total_successful_dribbles, 0)) as successful_dribbles,
        sum(coalesce(omp.total_progressive_run, 0)) as progressive_runs,
        sum(coalesce(omp.total_accelerations, 0)) as accelerations,

        -- defending and pressing
        sum(coalesce(omp.total_defensive_actions, 0)) as defensive_actions,
        sum(coalesce(omp.total_successful_defensive_actions, 0)) as successful_defensive_actions,
        sum(coalesce(omp.total_interceptions, 0)) as interceptions,
        sum(coalesce(omp.total_recoveries, 0)) as recoveries,
        sum(coalesce(omp.total_opponent_half_recoveries, 0)) as opponent_half_recoveries,
        sum(coalesce(omp.total_dangerous_opponent_half_recoveries, 0)) as dangerous_opponent_half_recoveries,
        sum(coalesce(omp.total_counterpressing_recoveries, 0)) as counterpressing_recoveries,
        sum(coalesce(omp.total_pressing_duels, 0)) as pressing_duels,
        sum(coalesce(omp.total_pressing_duels_won, 0)) as pressing_duels_won,
        sum(coalesce(omp.total_clearances, 0)) as clearances,
        sum(coalesce(omp.total_sliding_tackles, 0)) as sliding_tackles,
        sum(coalesce(omp.total_successful_sliding_tackles, 0)) as successful_sliding_tackles,

        -- duels and physical contests
        sum(coalesce(omp.total_duels, 0)) as duels,
        sum(coalesce(omp.total_duels_won, 0)) as duels_won,
        sum(coalesce(omp.total_defensive_duels, 0)) as defensive_duels,
        sum(coalesce(omp.total_defensive_duels_won, 0)) as defensive_duels_won,
        sum(coalesce(omp.total_offensive_duels, 0)) as offensive_duels,
        sum(coalesce(omp.total_offensive_duels_won, 0)) as offensive_duels_won,
        sum(coalesce(omp.total_aerial_duels, 0)) as aerial_duels,
        sum(coalesce(omp.total_aerial_duels_won, 0)) as aerial_duels_won,
        sum(coalesce(omp.total_field_aerial_duels, 0)) as field_aerial_duels,
        sum(coalesce(omp.total_field_aerial_duels_won, 0)) as field_aerial_duels_won,
        sum(coalesce(omp.total_loose_ball_duels, 0)) as loose_ball_duels,
        sum(coalesce(omp.total_loose_ball_duels_won, 0)) as loose_ball_duels_won,

        -- discipline and ball security
        sum(coalesce(omp.total_fouls, 0)) as fouls,
        sum(coalesce(omp.total_fouls_suffered, 0)) as fouls_suffered,
        sum(coalesce(omp.total_yellow_cards, 0)) as yellow_cards,
        sum(coalesce(omp.total_red_cards, 0)) as red_cards,
        sum(coalesce(omp.total_direct_red_cards, 0)) as direct_red_cards,
        sum(coalesce(omp.total_losses, 0)) as losses,
        sum(coalesce(omp.total_own_half_losses, 0)) as own_half_losses,
        sum(coalesce(omp.total_dangerous_own_half_losses, 0)) as dangerous_own_half_losses,
        sum(coalesce(omp.total_missed_balls, 0)) as missed_balls
    from outfield_match_player_base as omp
    group by
        omp.player_id,
        omp.season_id,
        omp.competition_id
),

/* ---------------------------------------------------------------------------
   5) Apply minimum exposure threshold for analytical stability
   --------------------------------------------------------------------------- */
player_season_filtered as (
    select
        psb.*
    from player_season_base as psb
    where psb.minutes >= 600
),

/* ---------------------------------------------------------------------------
   6) Engineer comparable player-season features
   --------------------------------------------------------------------------- */
player_season_features as (
    select
        psf.*,

        /* ---------------------------
           Per-90 intensity variables
           --------------------------- */
        90.0 * psf.goals / nullif(psf.minutes, 0) as goals_p90,
        90.0 * psf.assists / nullif(psf.minutes, 0) as assists_p90,
        90.0 * psf.shots / nullif(psf.minutes, 0) as shots_p90,
        90.0 * psf.shots_on_target / nullif(psf.minutes, 0) as shots_on_target_p90,
        90.0 * psf.xg_shot / nullif(psf.minutes, 0) as xg_p90,
        90.0 * psf.xg_assist / nullif(psf.minutes, 0) as xga_p90,
        90.0 * psf.touches_box / nullif(psf.minutes, 0) as touches_box_p90,

        90.0 * psf.passes / nullif(psf.minutes, 0) as passes_p90,
        90.0 * psf.progressive_passes / nullif(psf.minutes, 0) as progressive_passes_p90,
        90.0 * psf.passes_to_final_third / nullif(psf.minutes, 0) as passes_to_final_third_p90,
        90.0 * psf.key_passes / nullif(psf.minutes, 0) as key_passes_p90,
        90.0 * psf.shot_assists / nullif(psf.minutes, 0) as shot_assists_p90,

        90.0 * psf.dribbles / nullif(psf.minutes, 0) as dribbles_p90,
        90.0 * psf.progressive_runs / nullif(psf.minutes, 0) as progressive_runs_p90,
        90.0 * psf.accelerations / nullif(psf.minutes, 0) as accelerations_p90,

        90.0 * psf.defensive_actions / nullif(psf.minutes, 0) as defensive_actions_p90,
        90.0 * psf.interceptions / nullif(psf.minutes, 0) as interceptions_p90,
        90.0 * psf.recoveries / nullif(psf.minutes, 0) as recoveries_p90,
        90.0 * psf.opponent_half_recoveries / nullif(psf.minutes, 0) as opponent_half_recoveries_p90,
        90.0 * psf.dangerous_opponent_half_recoveries / nullif(psf.minutes, 0) as dangerous_opponent_half_recoveries_p90,
        90.0 * psf.counterpressing_recoveries / nullif(psf.minutes, 0) as counterpressing_recoveries_p90,
        90.0 * psf.pressing_duels / nullif(psf.minutes, 0) as pressing_duels_p90,
        90.0 * psf.clearances / nullif(psf.minutes, 0) as clearances_p90,
        90.0 * psf.sliding_tackles / nullif(psf.minutes, 0) as sliding_tackles_p90,

        90.0 * psf.duels / nullif(psf.minutes, 0) as duels_p90,
        90.0 * psf.defensive_duels / nullif(psf.minutes, 0) as defensive_duels_p90,
        90.0 * psf.offensive_duels / nullif(psf.minutes, 0) as offensive_duels_p90,
        90.0 * psf.aerial_duels / nullif(psf.minutes, 0) as aerial_duels_p90,
        90.0 * psf.field_aerial_duels / nullif(psf.minutes, 0) as field_aerial_duels_p90,
        90.0 * psf.loose_ball_duels / nullif(psf.minutes, 0) as loose_ball_duels_p90,

        /* ---------------------------
           Efficiency / success ratios
           --------------------------- */
        1.0 * psf.successful_passes / nullif(psf.passes, 0) as pass_accuracy,
        1.0 * psf.successful_smart_passes / nullif(psf.smart_passes, 0) as smart_pass_accuracy,
        1.0 * psf.successful_passes_to_final_third / nullif(psf.passes_to_final_third, 0) as final_third_pass_accuracy,
        1.0 * psf.successful_crosses / nullif(psf.crosses, 0) as cross_accuracy,
        1.0 * psf.successful_forward_passes / nullif(psf.forward_passes, 0) as forward_pass_accuracy,
        1.0 * psf.successful_through_passes / nullif(psf.through_passes, 0) as through_pass_accuracy,
        1.0 * psf.successful_key_passes / nullif(psf.key_passes, 0) as key_pass_accuracy,
        1.0 * psf.successful_vertical_passes / nullif(psf.vertical_passes, 0) as vertical_pass_accuracy,
        1.0 * psf.successful_long_passes / nullif(psf.long_passes, 0) as long_pass_accuracy,
        1.0 * psf.successful_progressive_passes / nullif(psf.progressive_passes, 0) as progressive_pass_accuracy,

        1.0 * psf.successful_dribbles / nullif(psf.dribbles, 0) as dribble_success_rate,
        1.0 * psf.successful_defensive_actions / nullif(psf.defensive_actions, 0) as defensive_action_success_rate,
        1.0 * psf.pressing_duels_won / nullif(psf.pressing_duels, 0) as pressing_duel_win_rate,
        1.0 * psf.duels_won / nullif(psf.duels, 0) as duel_win_rate,
        1.0 * psf.defensive_duels_won / nullif(psf.defensive_duels, 0) as defensive_duel_win_rate,
        1.0 * psf.offensive_duels_won / nullif(psf.offensive_duels, 0) as offensive_duel_win_rate,
        1.0 * psf.aerial_duels_won / nullif(psf.aerial_duels, 0) as aerial_duel_win_rate,
        1.0 * psf.field_aerial_duels_won / nullif(psf.field_aerial_duels, 0) as field_aerial_duel_win_rate,
        1.0 * psf.loose_ball_duels_won / nullif(psf.loose_ball_duels, 0) as loose_ball_duel_win_rate,

        1.0 * psf.shots_on_target / nullif(psf.shots, 0) as shot_on_target_rate,
        1.0 * psf.goals / nullif(psf.shots, 0) as goal_conversion_rate,
        1.0 * psf.xg_shot / nullif(psf.shots, 0) as xg_per_shot,

        /* ---------------------------
           Optional descriptive ratios
           --------------------------- */
        1.0 * psf.starts / nullif(psf.matches, 0) as start_rate,
        1.0 * psf.matches_substituted / nullif(psf.matches, 0) as substituted_rate,
        1.0 * psf.matches_coming_off / nullif(psf.matches, 0) as coming_off_rate,
        1.0 * psf.fouls_suffered / nullif(psf.fouls, 0) as fouls_suffered_per_foul_committed

    from player_season_filtered as psf
)

/* ---------------------------------------------------------------------------
   7) Final analytical output
   --------------------------------------------------------------------------- */
select
    psf.player_id,
    psf.player_name,
    psf.player_short_name,
    psf.player_name_unidecoded,
    psf.player_height,
    psf.player_weight,
    psf.player_birth_date,
    psf.player_birth_area_name,
    psf.player_passport_area_name,
    psf.role_name,
    psf.role_code_2,
    psf.role_code_3,
    psf.player_prefered_foot,

    psf.season_id,
    psf.season_name,
    psf.season_start_date,
    psf.season_end_date,
    psf.competition_id,
    psf.competition_name,
    psf.competition_area_name,
    psf.competition_format,
    psf.competition_type,
    psf.competition_category,

    psf.matches,
    psf.starts,
    psf.matches_substituted,
    psf.matches_coming_off,
    psf.minutes,
    psf.minutes_tagged,

    -- raw totals
    psf.goals,
    psf.assists,
    psf.shots,
    psf.shots_on_target,
    psf.xg_shot,
    psf.xg_assist,
    psf.passes,
    psf.progressive_passes,
    psf.dribbles,
    psf.defensive_actions,
    psf.duels,
    psf.aerial_duels,
    psf.pressing_duels,
    psf.recoveries,
    psf.opponent_half_recoveries,
    psf.touches_box,

    -- engineered features
    psf.goals_p90,
    psf.assists_p90,
    psf.shots_p90,
    psf.shots_on_target_p90,
    psf.xg_p90,
    psf.xga_p90,
    psf.touches_box_p90,

    psf.passes_p90,
    psf.progressive_passes_p90,
    psf.passes_to_final_third_p90,
    psf.key_passes_p90,
    psf.shot_assists_p90,

    psf.dribbles_p90,
    psf.progressive_runs_p90,
    psf.accelerations_p90,

    psf.defensive_actions_p90,
    psf.interceptions_p90,
    psf.recoveries_p90,
    psf.opponent_half_recoveries_p90,
    psf.dangerous_opponent_half_recoveries_p90,
    psf.counterpressing_recoveries_p90,
    psf.pressing_duels_p90,
    psf.clearances_p90,
    psf.sliding_tackles_p90,

    psf.duels_p90,
    psf.defensive_duels_p90,
    psf.offensive_duels_p90,
    psf.aerial_duels_p90,
    psf.field_aerial_duels_p90,
    psf.loose_ball_duels_p90,

    psf.pass_accuracy,
    psf.smart_pass_accuracy,
    psf.final_third_pass_accuracy,
    psf.cross_accuracy,
    psf.forward_pass_accuracy,
    psf.through_pass_accuracy,
    psf.key_pass_accuracy,
    psf.vertical_pass_accuracy,
    psf.long_pass_accuracy,
    psf.progressive_pass_accuracy,

    psf.dribble_success_rate,
    psf.defensive_action_success_rate,
    psf.pressing_duel_win_rate,
    psf.duel_win_rate,
    psf.defensive_duel_win_rate,
    psf.offensive_duel_win_rate,
    psf.aerial_duel_win_rate,
    psf.field_aerial_duel_win_rate,
    psf.loose_ball_duel_win_rate,

    psf.shot_on_target_rate,
    psf.goal_conversion_rate,
    psf.xg_per_shot,

    psf.start_rate,
    psf.substituted_rate,
    psf.coming_off_rate,
    psf.fouls_suffered_per_foul_committed

from player_season_features as psf
order by
    psf.season_start_date,
    psf.competition_id,
    psf.player_name;
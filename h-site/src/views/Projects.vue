<template>
  <v-container v-if="projects.length > 0">
    <v-row v-for="(project, idx) in projects" :key="idx">
      <v-col>{{ project.name }}</v-col>
    </v-row>
  </v-container>
</template>

<script>
import { mapActions, mapGetters } from 'vuex'

export default {
  name: "Projects",
  methods: {
    ...mapActions(['getProjects'])
  },
  computed: {
    ...mapGetters(['projects'])
  },
  created() {
    this.getProjects().then(console.log("Projects fetched"))
  }
};
</script>

<style></style>

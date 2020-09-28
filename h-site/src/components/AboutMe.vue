<template>
  <v-container fluid>
    <v-row>
      <v-col v-for="(obj, idx) in Om" :key="idx" class="text-start" md="4" sm="8">
        <h3>{{ obj.title }}</h3>
        <v-list v-for="(element, index) in obj.desc" :key="index">{{ element }}</v-list>
      </v-col>
    </v-row>
  </v-container>
</template>

<script>
export default {
  data: function() {
    return {
      Om: [
        {
          id: 1,
          title: "About me",
          desc: [
            "After finishing my MSc in Astrophysics at UiO I started to play around with HTML and css while looking for a job.",
            "Quite early I knew that this was something I wanted to master and work with, so I reached out to a close friend in the business and asked what I should learn.",
            "He recommended learning a frontend framework, and as a fan he recommended Vuejs.",
            "Thus, last summer I started to learn Vuejs and decided to take an extra year at university taking programming classes to get a more fundamental understanding of programming.",
            "A little info about this site: It's built with Vuejs, using vuex, vue-router and vuetify. In the backend I have Firestore. You need to sign in to use the todo manager and profile view, and this is a feature reserved for myself at the time being. But feel free to have a look around!"
          ]
        },
        {
          id: 2,
          title: "Degrees",
          desc: [
            "MSc, Astronomy (Cosmology) at UiO",
            "BSc, Physics (Astronomy) at UiO"
          ]
        },
        {
          id: 3,
          title: "The Last Year",
          desc: [
            "Functional Programming (scheme) - UiO",
            "OOP (Java) - UiO",
            "Databases and data modelling (relational databases, SQL) - UiO",
            "Databasesystems (SQL, MongoDB++) - UiO",
            "High Performance Computing (C) - UiO",
            "Introductory Computertechnology - UiO",
            "Vuejs + libraries",
            "Firebase Firestore"
          ]
        }
      ]
    };
  }
};
</script>

<style scoped>
h3 {
  font-size: 130%;
}
</style>
